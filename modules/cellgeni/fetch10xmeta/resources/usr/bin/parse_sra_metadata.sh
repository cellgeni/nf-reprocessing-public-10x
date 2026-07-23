#!/bin/bash

SERIES=$1

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

for i in `cat $SERIES.run.list`
do
  TYPE="SRA"  ## we always default to SRA. This could cause problems for very fresh datasets. 
  LOC=""

  SPECIES=`grep -w $i $SERIES.sra.tsv | cut -f29`
  SPECIES=${SPECIES:-UNKNOWN}
  SRA=`grep -w $i $SERIES.sra.tsv | cut -f10`

  SUCCESS=0

  # Try getting BAM file from SDL api
  SDLBAM=`sdl_retrieve "$i" | jq -r '
          .result[].files[] |
          select(.name | contains("bam")) |
          .locations[] |
          select((.rehydrationRequired // false) == false and (.payRequired // false) == false) |
          .link
        '`
  sleep 0.3
  if [[ $SDLBAM != "" ]]
  then
    TYPE="BAM"
    LOC=$SDLBAM
    echo $SDLBAM >> $SERIES.urls.list
    >&2 echo "Sample $i is available via SRA as an original submitter's BAM file: $LOC"
    SUCCESS=1
  fi

  # Try getting SRA file from SDL api
  if [[ $SUCCESS -eq 0 ]]
  then
    SDLSRA=`sdl_retrieve "$i" | jq -r '
              [
                .result[]
                .files[]
                | select(.type == "sra")
                | .locations[]
                | select((.payRequired // false) == false)
              ]
              | map(.link)
              | first'`
    sleep 0.3
    if [[ $SDLSRA != "" ]]
    then
      TYPE="SRA"
      LOC=$SDLSRA
      echo $SDLSRA >> $SERIES.urls.list
      >&2 echo "Sample $i is available via NCBI/Amazon as an SRA archive: $LOC"
      SUCCESS=1
    elif [[ $SRA != "" ]]
    then
      LOC=$SRA
      echo $SRA >> $SERIES.urls.list
      >&2 echo "Sample $i is available via SRA as an SRA archive: $LOC"
      SUCCESS=1
    else
      SRA=`srapath $i`
      LOC=$SRA
      >&2 echo "WARNING: No ENA ftp URL found for sample $i, using 'srapath' to get the (open) Amazon link to SRA archive.."
      echo $SRA >> $SERIES.urls.list
      >&2 echo "Sample $i is available via NCBI/Amazon as an SRA archive: $LOC"
    fi
  fi

  echo -e "$i\t$SPECIES\t$LOC\t$TYPE"
done

