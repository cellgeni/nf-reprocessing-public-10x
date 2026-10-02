#!/usr/bin/env python3
"""
Send a batch-debug notification by email and, optionally, Slack.

Email needs no credentials: the farm's local Postfix relay accepts mail from
any node with no auth, so this connects to smtplib.SMTP("localhost") exactly
like track_reprocessing/src/tracking/cli.py's emit_reports() does. There is
nothing to configure and nothing this script can get wrong about auth.

Slack is optional and does need a credential — SLACK_BOT_TOKEN and CHANNEL_ID,
the same two variables reprocessing_slack_bot/app.py reads. This script never
invents or guesses them: with no --slack flag, or with the flag set but no
token resolvable, it skips Slack and says so on stderr, and the email still
sends. It never fails the whole run just because Slack isn't wired up.

Slack's one library dependency, slack_sdk, is satisfied through `uvx` rather
than an installed package: the farm python has user site-packages DISABLED
(site.ENABLE_USER_SITE is False), so `pip install --user slack_sdk` writes a
directory that nothing will ever import, which is how this silently skipped
Slack for the batch 8 notification on 2026-09-16. uv and uvx are on PATH at
/software/cellgen/cellgeni/uv/. An in-process import is still tried first, so
an environment that does have slack_sdk pays nothing for this.
"""
import argparse
import os
import re
import shutil
import smtplib
import subprocess
import sys
from email.message import EmailMessage
from pathlib import Path

DEFAULT_FROM = "noreply-reprocessing@cellgeni-su"

# Everyone who should get a batch notification. This list is the single place to
# add or drop a recipient: --team expands to it, and the commands in SKILL.md §8
# and references/notify.md §3 use --team rather than spelling addresses out, so a
# name added here reaches every future notification without editing any prose.
# ap41 was in none of those invocations until 2026-09-22 and therefore silently
# received nothing for batches 6-9 — the mail relay was never the problem.
TEAM = ["ab76@sanger.ac.uk", "ap41@sanger.ac.uk"]


def resolve_recipients(to: list[str], team: bool) -> list[str]:
    """Flatten comma- or space-separated --to values, add TEAM if asked, dedupe.

    Accepting "a@x,b@y" as a single token matters because notify_run_done.sh's
    --to takes one shell argument and passes it straight through.
    """
    out = []
    for item in ([*to, *TEAM] if team else to):
        for addr in re.split(r"[,\s]+", item):
            if addr and addr not in out:
                out.append(addr)
    return out


def load_env_file(path: str) -> dict:
    """Minimal KEY=VALUE .env reader — avoids a hard dependency on python-dotenv."""
    env = {}
    p = Path(path)
    if not p.exists():
        return env
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        env[key.strip()] = value.strip().strip('"').strip("'")
    return env


def send_email(to: list[str], subject: str, body: str, attach: list[str], sender: str) -> None:
    msg = EmailMessage()
    msg.set_content(body)
    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = ", ".join(to)
    for path in attach:
        p = Path(path)
        if not p.is_file():
            print(f"warn: --attach {path} not found, skipping", file=sys.stderr)
            continue
        msg.add_attachment(p.read_bytes(), maintype="application", subtype="octet-stream", filename=p.name)
    with smtplib.SMTP("localhost") as server:
        server.send_message(msg)
    print(f"emailed {', '.join(to)}: {subject}", file=sys.stderr)


# Run inside the ephemeral uvx environment. Token and channel arrive through the
# child's environment, never argv: /proc/<pid>/environ is owner-readable while
# /proc/<pid>/cmdline is world-readable, and this runs on shared login nodes.
_SLACK_SNIPPET = r"""
import os, sys
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
try:
    WebClient(token=os.environ["SLACK_BOT_TOKEN"]).chat_postMessage(
        channel=os.environ["CHANNEL_ID"], text=sys.stdin.read())
except SlackApiError as e:
    sys.stderr.write(e.response["error"]); sys.exit(2)
"""


def _post_via_uvx(text: str, token: str, channel: str) -> bool:
    """slack_sdk is not importable here — borrow it for one call via uvx."""
    uvx = shutil.which("uvx") or shutil.which("uv")
    if not uvx:
        print("slack: slack_sdk not importable and no uvx/uv on PATH — skipped "
              "(module load cellgen/uv, or pip install slack_sdk into a venv)", file=sys.stderr)
        return False
    cmd = ([uvx, "--with", "slack_sdk", "python", "-c", _SLACK_SNIPPET]
           if Path(uvx).name == "uvx" else
           [uvx, "run", "--no-project", "--with", "slack_sdk", "python", "-c", _SLACK_SNIPPET])
    env = {**os.environ, "SLACK_BOT_TOKEN": token, "CHANNEL_ID": channel}
    try:
        proc = subprocess.run(cmd, input=text, text=True, env=env,
                              capture_output=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"slack: {Path(uvx).name} could not run slack_sdk, skipped — {e}", file=sys.stderr)
        return False
    if proc.returncode == 0:
        print(f"slack: posted to channel {channel} (slack_sdk via {Path(uvx).name})", file=sys.stderr)
        return True
    detail = (proc.stderr or proc.stdout or "").strip().splitlines()
    print(f"slack: skipped — {detail[-1] if detail else f'exit {proc.returncode}'}", file=sys.stderr)
    return False


def send_slack(text: str, token: str | None, channel: str | None) -> bool:
    """Returns True if actually sent. Never raises — a missing token or a Slack
    API error is reported on stderr and treated as 'skipped', not fatal."""
    if not token or not channel:
        print("slack: SLACK_BOT_TOKEN/CHANNEL_ID not set — skipped (email still sent)", file=sys.stderr)
        return False
    try:
        from slack_sdk import WebClient
        from slack_sdk.errors import SlackApiError
    except ImportError:
        return _post_via_uvx(text, token, channel)
    try:
        client = WebClient(token=token)
        client.chat_postMessage(channel=channel, text=text)
        print(f"slack: posted to channel {channel}", file=sys.stderr)
        return True
    except SlackApiError as e:
        print(f"slack: API error, skipped — {e.response['error']}", file=sys.stderr)
        return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--to", nargs="+", default=[],
                    help="Recipient email address(es); a comma-separated list in a single "
                         "argument works too. Omit for a Slack-only message (--slack with "
                         "--slack-text/--slack-text-file), which is what notify_run_done.sh does")
    ap.add_argument("--team", action="store_true",
                    help="Add the standard batch-notification recipients (" + ", ".join(TEAM) +
                         ") to --to. This is the intended way to send a batch notification — "
                         "spelling a single address out by hand is how a recipient goes missing")
    ap.add_argument("--subject", default=None, help="Required when --to is given")
    ap.add_argument("--body", default=None, help="Message body text (mutually exclusive with --body-file)")
    ap.add_argument("--body-file", default=None)
    ap.add_argument("--from", dest="sender", default=DEFAULT_FROM)
    ap.add_argument("--attach", nargs="*", default=[])
    ap.add_argument("--slack", action="store_true", help="Also post to Slack (mrkdwn is a different language from the email body — see notify.md's two templates)")
    ap.add_argument("--slack-text", default=None, help="Slack mrkdwn text, inline. Falls back to the email body if neither this nor --slack-text-file is given")
    ap.add_argument("--slack-text-file", default=None, help="Slack mrkdwn text, from a file (preferred over --slack-text for anything multi-line)")
    ap.add_argument("--slack-env-file", default=None,
                     help="Optional .env-style file to read SLACK_BOT_TOKEN/CHANNEL_ID from, "
                          "if they are not already in the environment")
    ap.add_argument("--dry-run", action="store_true", help="Print what would be sent, send nothing")
    args = ap.parse_args()
    recipients = resolve_recipients(args.to, args.team)

    if not recipients and not args.slack:
        ap.error("nothing to send: give --to/--team for email, --slack for Slack, or both")
    if recipients and not args.subject:
        ap.error("--subject is required with --to/--team")
    if args.body and args.body_file:
        ap.error("--body and --body-file are mutually exclusive")

    if args.slack_text_file and args.slack_text:
        ap.error("--slack-text and --slack-text-file are mutually exclusive")
    slack_text = args.slack_text
    if args.slack_text_file:
        slack_text = Path(args.slack_text_file).read_text()

    body = args.body if args.body is not None else (
        Path(args.body_file).read_text() if args.body_file else None)
    if recipients and body is None:
        ap.error("exactly one of --body or --body-file is required with --to/--team")
    if args.slack and slack_text is None and body is None:
        ap.error("--slack with no --to needs --slack-text or --slack-text-file")
    slack_text = slack_text or body  # only a same-text fallback if neither was given at all

    if args.dry_run:
        if recipients:
            print(f"--- would email {', '.join(recipients)} ---\nSubject: {args.subject}\n\n{body}\n---")
        if args.slack:
            print(f"--- would post to Slack ---\n{slack_text}\n---")
        return 0

    if recipients:
        send_email(recipients, args.subject, body, args.attach, args.sender)

    if args.slack:
        env = {}
        if args.slack_env_file:
            env = load_env_file(args.slack_env_file)
        token = os.environ.get("SLACK_BOT_TOKEN") or env.get("SLACK_BOT_TOKEN")
        channel = os.environ.get("CHANNEL_ID") or env.get("CHANNEL_ID")
        send_slack(slack_text, token, channel)

    return 0


if __name__ == "__main__":
    sys.exit(main())
