#!/bin/bash
#
# Work out, for every run of a series, which species it is and where its reads
# can be downloaded from. Replaces parse_ena_metadata.sh and
# parse_sra_metadata.sh, each of which read a single table and were called
# either/or. Reading both together fixes two things:
#
#   * Species is taken from ENA and SRA together. A run is often listed in only
#     one of the tables — in GSE166796 every run absent from ena.tsv had a
#     perfectly good species in sra.tsv — and reading one table alone reported
#     those runs as 'UNKNOWN'. Downstream that made the runs of a single sample
#     disagree about their species, which split the sample in two.
#   * The NCBI SDL API is queried at most once per run. The BAM lookup and the
#     SRA lookup previously issued a separate request for the same accession.
#
# Download location is picked by this order of preference:
#
#    1. ArrayExpress submitter fastq  (sdrf.txt, E-MTAB series only, and only
#                                      files the study actually registers)
#    2. ENA paired-end fastq          (ena.tsv col 11, strict _1/_2 naming)
#    3. ENA submitter fastq           (ena.tsv col 12, when no BAM is offered)
#    4. ENA submitter BAM             (ena.tsv col 12)
#    5. SDL BAM                       (the single API call)
#    6. ENA SRA archive               (ena.tsv col 13)
#    7. SDL SRA archive               (reuses the cached API response)
#    8. SRA archive                   (sra.tsv col 10)
#    9. srapath
#
# This is the two original orderings preserved as-is: the ENA parser preferred
# its own col 13 over SDL, the SRA parser preferred SDL over its col 10.
# Steps 1-4 are decided from the tables alone, so runs that get that far never
# reach the network here at all.
#
# Usage: parse_metadata.sh <series_id>
#   reads  <series>.run.list and whichever of <series>.ena.tsv, <series>.sra.tsv,
#          <series>.sdrf.txt and <series>.aefiles.list exist
#   writes <series>.urls.list
#   prints run <TAB> species <TAB> location <TAB> type

set -uo pipefail

SERIES=$1

if [[ ! -s $SERIES.run.list ]]
then
  >&2 echo "ERROR: No run list '$SERIES.run.list' found!"
  exit 1
fi

if [[ ! -s $SERIES.ena.tsv && ! -s $SERIES.sra.tsv ]]
then
  >&2 echo "ERROR: Neither '$SERIES.ena.tsv' nor '$SERIES.sra.tsv' found!"
  exit 1
fi

if [[ -f $SERIES.urls.list ]]
then
  >&2 echo "WARNING: File '$SERIES.urls.list' exists! This should not happen; overwriting the file.."
  rm $SERIES.urls.list
fi

# Query the NCBI SDL API for one accession, printing the raw JSON body on stdout.
# Retries on transient network/HTTP failures or an empty response body; returns
# non-zero (and empty output) only after all attempts are exhausted, so callers
# fall back exactly as they would on a genuinely empty result.
sdl_retrieve() {
  local acc=$1
  local url="https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=${acc}&accept-alternate-locations=yes"
  local tries=${SDL_MAX_TRIES:-4}
  local attempt body
  for attempt in $(seq 1 "$tries")
  do
    if body=$(curl -sS --fail --retry 3 --retry-delay 2 --retry-connrefused --max-time 60 "$url" 2>/dev/null) && [[ -n $body ]]
    then
      printf '%s' "$body"
      return 0
    fi
    >&2 echo "WARNING: SDL query for $acc failed (attempt $attempt/$tries); retrying in $((attempt * 2))s.."
    sleep $((attempt * 2))
  done
  >&2 echo "WARNING: SDL query for $acc failed after $tries attempts; continuing with empty response"
  return 1
}

# Fetch the SDL response for one run at most once and leave it in RUN_SDL_JSON,
# which the BAM and the SRA lookup then both read.
#
# This deliberately returns nothing and must never be used as a pipeline stage:
# bash runs those in a subshell, so the cached body would be thrown away the
# moment it returned and every lookup would issue its own request — precisely
# the duplication this replaces.
RUN_SDL_JSON=""
RUN_SDL_TRIED=0
sdl_fetch() {
  if [[ $RUN_SDL_TRIED -eq 0 ]]
  then
    RUN_SDL_TRIED=1
    RUN_SDL_JSON=$(sdl_retrieve "$1") || RUN_SDL_JSON=""
    sleep 0.3
  fi
}

# Run a jq filter over the cached SDL body, quietly yielding nothing when the
# API gave us no usable answer.
sdl_query() {
  [[ -n $RUN_SDL_JSON ]] || return 0
  printf '%s' "$RUN_SDL_JSON" | jq -r "$1" 2>/dev/null
}

# Join the tables and decide everything that does not need the network. Emits
# one line per run of <series>.run.list:
#
#   run <TAB> species <TAB> location <TAB> type <TAB> ena_sra <TAB> sra_url
#
# '-' stands for an absent value, so that no field is ever empty — bash treats
# runs of tabs as a single separator, which would silently shift the columns.
join_tables() {
  awk -F'\t' -v OFS='\t' \
      -v RUNS="$SERIES.run.list"       -v ENA="$SERIES.ena.tsv" \
      -v SRA="$SERIES.sra.tsv"         -v SDRF="$SERIES.sdrf.txt" \
      -v AEFILES="$SERIES.aefiles.list" \
      -v AEBASE="https://www.ebi.ac.uk/biostudies/files/$SERIES/" \
      -v AEMIRROR="/pub/databases/(microarray|arrayexpress)/data/experiment/" '
    function dash(s) { return (s == "" ? "-" : s) }

    # Add the fastq URLs of one SDRF row to a run. A run can span several rows
    # — one per lane, or one per read file — so these accumulate rather than
    # overwrite, and ae_seen keeps a URL repeated across rows from being listed
    # (and downloaded) twice.
    function ae_add(r, list,   k, u, i) {
      if (list == "") return
      k = split(list, u, ";")
      for (i = 1; i <= k; i++) {
        if ((r SUBSEP u[i]) in ae_seen) continue
        ae_seen[r SUBSEP u[i]] = 1
        ae_fq[r] = (ae_fq[r] == "" ? u[i] : ae_fq[r] ";" u[i])
      }
    }

    # run list first, so the ArrayExpress pass below knows which runs to look for
    FILENAME == RUNS { if ($1 != "") { order[++n] = $1; runs[$1] = 1 } next }

    # ena.tsv: 1=run 10=species 11=fastq_ftp 12=submitted_ftp 13=sra_ftp
    FILENAME == ENA { ena_sp[$1]=$10; ena_fq[$1]=$11; ena_sub[$1]=$12; ena_sra[$1]=$13; next }

    # sra.tsv: 1=run 10=archive url 29=species
    FILENAME == SRA { sra_sp[$1]=$29; sra_url[$1]=$10; next }

    # the files the study registers, one base name per line, absent when
    # collect_metadata.sh could not reach the BioStudies API
    FILENAME == AEFILES { if ($1 != "") { aefiles[$1] = 1; have_aefiles = 1 } next }

    # sdrf.txt: submitter fastq URLs sit in arbitrary columns of the run row.
    #
    # A URI under the pre-BioStudies mirror at
    # ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/ is checked against
    # the files the study registers before it is believed, because that spelling
    # means two different things. E-MTAB-8060 and E-MTAB-9221 look identical
    # here — both point every Comment[FASTQ_URI] at the mirror, and both also
    # declare a Comment[BAM_URI] for a registered ENA submission — yet:
    #
    #   * E-MTAB-9221 registers all 40 of those fastq files. The URI is a stale
    #     spelling of a real study file, so it is kept and the URL refreshed.
    #   * E-MTAB-8060 registers none of its 36. The study holds only the idf and
    #     the sdrf; the reads exist solely on the unmaintained mirror, and the
    #     real submission is the ENA BAM. Taking those URIs at face value
    #     discarded the BAM for all 15 runs, silently, because rule 1 never
    #     looks at what it is outranking. Dropping them leaves ae_fq empty, and
    #     the run falls through to the ENA table and picks the BAM up.
    #
    # Anything not on the mirror is left exactly as the SDRF wrote it, and so is
    # everything when the file list is missing: an unreachable API should not be
    # able to reroute a whole series to its BAMs.
    FILENAME == SDRF {
      fq = ""; unreg = ""
      for (i = 1; i <= NF; i++) {
        if ($i !~ /ftp:\/\/.*\.f.*q/) continue
        if (!have_aefiles || $i !~ AEMIRROR) { fq = (fq == "" ? $i : fq ";" $i); continue }
        base = $i
        sub(/.*\//, "", base)
        if (base in aefiles)
          fq = (fq == "" ? AEBASE base : fq ";" AEBASE base)
        else
          unreg = (unreg == "" ? $i : unreg ";" $i)
      }
      if (fq == "" && unreg == "") next
      for (r in runs) {
        if ($0 !~ ("(^|[^A-Za-z0-9])" r "([^A-Za-z0-9]|$)")) continue
        ae_add(r, fq)
        if (unreg != "") ae_unreg[r] += split(unreg, discarded, ";")
      }
      next
    }

    END {
      for (j = 1; j <= n; j++) {
        r = order[j]

        # species: whichever table has a real answer, ENA first
        sp = ena_sp[r]
        if (sp == "" || sp == "UNKNOWN") sp = sra_sp[r]
        if (sp == "") sp = "UNKNOWN"

        # ENA insists on _1/_2 naming for the reads it serves itself
        enagz = (ena_fq[r] ~ /_1\.fastq\.gz/ && ena_fq[r] ~ /_2\.fastq\.gz/) ? ena_fq[r] : ""

        # submitter uploads: anything fastq-ish, and BAMs minus their indexes,
        # which are very often uploaded alongside
        sub_ftp = ena_sub[r]
        orifq = (sub_ftp ~ /f.*q/) ? sub_ftp : ""
        oribam = ""
        m = split(sub_ftp, parts, ";")
        for (i = 1; i <= m; i++)
          if (parts[i] ~ /\.bam/ && parts[i] !~ /\.bai/)
            oribam = (oribam == "" ? parts[i] : oribam ";" parts[i])

        # say so out loud: the whole reason E-MTAB-8060 went unnoticed is that
        # the substitution left no trace anywhere
        if (ae_unreg[r] > 0)
          print "WARNING: run " r ": ignored " ae_unreg[r] \
                " SDRF fastq URI(s) naming files the study does not register" \
                (ae_fq[r] == "" ? "; using the ENA table instead" : "") \
                > "/dev/stderr"

        loc = ""; type = ""
        if (ae_fq[r] != "")                 { type = "ORIFQ"; loc = ae_fq[r] }
        else if (enagz != "")               { type = "ENAFQ"; loc = enagz }
        else if (orifq != "" && oribam == "") { type = "ORIFQ"; loc = orifq }
        else if (oribam != "")              { type = "BAM";   loc = oribam }

        # both kept separate: they sit either side of the SDL SRA lookup
        print r, sp, dash(loc), dash(type), dash(ena_sra[r]), dash(sra_url[r])
      }
    }
  ' "$SERIES.run.list" \
    $( [[ -s $SERIES.aefiles.list ]] && echo "$SERIES.aefiles.list" ) \
    $( [[ -s $SERIES.ena.tsv  ]] && echo "$SERIES.ena.tsv" ) \
    $( [[ -s $SERIES.sra.tsv  ]] && echo "$SERIES.sra.tsv" ) \
    $( [[ -s $SERIES.sdrf.txt ]] && echo "$SERIES.sdrf.txt" )
}

while IFS=$'\t' read -r RUN SPECIES LOC TYPE ENA_SRA SRA_URL
do
  [[ $LOC     == "-" ]] && LOC=""
  [[ $TYPE    == "-" ]] && TYPE=""
  [[ $ENA_SRA == "-" ]] && ENA_SRA=""
  [[ $SRA_URL == "-" ]] && SRA_URL=""

  RUN_SDL_JSON=""
  RUN_SDL_TRIED=0

  # 1-4: already settled from the tables
  if [[ -n $LOC ]]
  then
    echo "$LOC" | tr ';' '\n' >> $SERIES.urls.list
    >&2 echo "Run $RUN is available as $TYPE: $LOC"
  else
    # 5: a 10x BAM keeps the original reads, so it beats any SRA archive
    sdl_fetch "$RUN"
    SDLBAM=$(sdl_query '
            .result[].files[] |
            select(.name | contains("bam")) |
            .locations[] |
            select((.rehydrationRequired // false) == false and (.payRequired // false) == false) |
            .link
          ')
    if [[ -n $SDLBAM ]]
    then
      TYPE="BAM"
      LOC=$SDLBAM
      echo "$SDLBAM" | tr ';' '\n' >> $SERIES.urls.list
      >&2 echo "Run $RUN is available via NCBI/Amazon as a BAM file: $LOC"
    # 6: ENA's own SRA mirror
    elif [[ -n $ENA_SRA ]]
    then
      TYPE="SRA"
      LOC=$ENA_SRA
      echo "$ENA_SRA" >> $SERIES.urls.list
      >&2 echo "Run $RUN is available via ENA as an SRA archive: $LOC"
    else
      # 7: same SDL response as the BAM lookup above, no second request.
      # 'first' on an empty array yields the string "null", which the two
      # scripts this replaces would have happily used as a URL.
      SDLSRA=$(sdl_query '
                [
                  .result[]
                  .files[]
                  | select(.type == "sra")
                  | .locations[]
                  | select((.payRequired // false) == false)
                ]
                | map(.link)
                | first')
      if [[ -n $SDLSRA && $SDLSRA != "null" ]]
      then
        TYPE="SRA"
        LOC=$SDLSRA
        echo "$SDLSRA" >> $SERIES.urls.list
        >&2 echo "Run $RUN is available via NCBI/Amazon as an SRA archive: $LOC"
      # 8: the archive URL recorded in the SRA table
      elif [[ -n $SRA_URL ]]
      then
        TYPE="SRA"
        LOC=$SRA_URL
        echo "$SRA_URL" >> $SERIES.urls.list
        >&2 echo "Run $RUN is available via SRA as an SRA archive: $LOC"
      else
        # 9: neither table nor SDL knew; ask the toolkit directly
        TYPE="SRA"
        LOC=`srapath $RUN`
        >&2 echo "WARNING: No table or SDL URL found for run $RUN, using 'srapath' to get the (open) Amazon link to SRA archive.."
        echo "$LOC" >> $SERIES.urls.list
        >&2 echo "Run $RUN is available via NCBI/Amazon as an SRA archive: $LOC"
      fi
    fi
  fi

  echo -e "$RUN\t$SPECIES\t$LOC\t$TYPE"
done < <(join_tables)
