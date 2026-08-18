#!/usr/bin/env python3
"""Compile samtools flagstat reports into a single wide TSV.

Reads samtools' DEFAULT flagstat output format (not `-O tsv`). The mapping rule emits
the default format because MultiQC detects flagstat files by content, matching the
string "in total (QC-passed reads + QC-failed reads)" — the "in " prefix only appears in
the default format, so a tsv-format file would be silently ignored by the report. Since
the default (streaming) mapping path gets only one pass over the alignments, that same
file has to serve both MultiQC and this summary.

Each default-format line looks like:

    123456 + 0 in total (QC-passed reads + QC-failed reads)
    120000 + 0 mapped (97.24% : N/A)

i.e. "<QC-passed> + <QC-failed> <label>", where the label may carry a trailing
percentage in parentheses. Percentages are recomputed here rather than scraped from that
text, so that zero denominators yield 'NA' instead of a nonsense value.
"""

import argparse
import csv
import os
import re
import sys


STATUS_OK = 'ok'
STATUS_EMPTY = 'empty_input'

# "<passed> + <failed> <label>", label optionally followed by " (nn.nn% : N/A)"
FLAGSTAT_LINE = re.compile(r'^(\d+)\s*\+\s*(\d+)\s+(.+?)\s*$')

# samtools appends its percentage annotation as "(97.24% : N/A)", and "(N/A : N/A)" when
# the denominator is zero — so keying on '%' would miss the zero case. The ' : ' separator
# is the reliable signature. Crucially this leaves "(mapQ>=5)" alone, which carries meaning
# and distinguishes that metric from its non-mapQ counterpart.
TRAILING_PCT = re.compile(r'\s*\([^)]*:[^)]*\)\s*$')

# Metric labels as samtools writes them, mapped to output column names. Column order
# below defines the TSV column order.
COUNT_COLUMNS = [
    ('total',                                        'total'),
    ('primary',                                      'primary'),
    ('secondary',                                    'secondary'),
    ('supplementary',                                'supplementary'),
    ('duplicates',                                   'duplicates'),
    ('mapped',                                       'mapped'),
    ('primary mapped',                               'primary_mapped'),
    ('paired in sequencing',                         'paired_in_sequencing'),
    ('read1',                                        'read1'),
    ('read2',                                        'read2'),
    ('properly paired',                              'properly_paired'),
    ('with itself and mate mapped',                  'with_itself_and_mate_mapped'),
    ('singletons',                                   'singletons'),
    ('with mate mapped to a different chr',          'mate_mapped_diff_chr'),
    ('with mate mapped to a different chr (mapQ>=5)', 'mate_mapped_diff_chr_mapq5'),
]

# (percentage column, numerator column, denominator column)
PCT_COLUMNS = [
    ('mapped_pct',           'mapped',          'total'),
    ('primary_mapped_pct',   'primary_mapped',  'primary'),
    ('properly_paired_pct',  'properly_paired', 'paired_in_sequencing'),
    ('singletons_pct',       'singletons',      'paired_in_sequencing'),
]

# Final column order, matching the standalone mapNstat.sh summary plus a status column
HEADER = [
    'sample', 'total', 'primary', 'secondary', 'supplementary', 'duplicates',
    'mapped', 'mapped_pct', 'primary_mapped', 'primary_mapped_pct',
    'paired_in_sequencing', 'read1', 'read2',
    'properly_paired', 'properly_paired_pct',
    'with_itself_and_mate_mapped', 'singletons', 'singletons_pct',
    'mate_mapped_diff_chr', 'mate_mapped_diff_chr_mapq5',
    'status',
]


def normalise_label(label):
    """Reduce a raw flagstat line label to the key used in COUNT_COLUMNS.

    Strips samtools' percentage annotation, then folds the total line — written as
    "in total (QC-passed reads + QC-failed reads)", whose parenthesised part contains no
    colon and so survives the strip — down to plain "total".
    """
    label = TRAILING_PCT.sub('', label).strip()
    if label.startswith('in total'):
        return 'total'
    return label


def pct(numerator, denominator):
    """Percentage to 2dp, or 'NA' when the denominator is zero or unavailable.

    An unmapped or empty sample legitimately has a zero denominator, so this must not
    raise — it mirrors the pct() helper in the standalone mapNstat.sh script.
    """
    try:
        denominator = float(denominator)
        numerator = float(numerator)
    except (TypeError, ValueError):
        return 'NA'
    if denominator <= 0:
        return 'NA'
    return f"{(numerator / denominator) * 100:.2f}"


def empty_sample_stats():
    """Row for a sample whose flagstat file is a zero-byte placeholder.

    rule reference_mapping touches empty outputs when fastp produced no usable reads.
    Such a sample must still appear in the summary as zero rather than vanishing.
    """
    stats = {col: 0 for _, col in COUNT_COLUMNS}
    stats.update({col: 'NA' for col, _, _ in PCT_COLUMNS})
    stats['status'] = STATUS_EMPTY
    return stats


def parse_flagstat(path):
    """Parse one default-format flagstat file.

    Returns a stats dict, or None if the file is non-empty but unparseable.
    A zero-byte file is not an error: it is the placeholder for a sample with no reads.
    """
    try:
        if os.path.getsize(path) == 0:
            print(f"Note: {path} is empty (sample had no usable reads); "
                  f"recording zero counts", file=sys.stderr)
            return empty_sample_stats()

        found = {}
        with open(path, 'r') as fh:
            for line in fh:
                match = FLAGSTAT_LINE.match(line)
                if not match:
                    continue
                qc_passed, _qc_failed, label = match.groups()
                found[normalise_label(label)] = int(qc_passed)

        if not found:
            print(f"ERROR: no flagstat metrics found in {path}", file=sys.stderr)
            return None

        stats = {}
        missing = []
        for label, column in COUNT_COLUMNS:
            if label in found:
                stats[column] = found[label]
            else:
                stats[column] = 'NA'
                missing.append(label)

        # 'total' and 'primary' are present in every samtools version that emits this
        # format; anything else missing is a version difference, not a broken file.
        for required in ('total', 'primary'):
            if stats[required] == 'NA':
                print(f"ERROR: {path} is missing the '{required}' metric", file=sys.stderr)
                return None

        if missing:
            print(f"Note: {path} has no {', '.join(missing)} line(s); "
                  f"recorded as NA", file=sys.stderr)

        for column, numerator, denominator in PCT_COLUMNS:
            stats[column] = pct(stats[numerator], stats[denominator])

        stats['status'] = STATUS_OK
        return stats

    except (OSError, ValueError) as e:
        print(f"ERROR: could not parse {path}: {e}", file=sys.stderr)
        return None


def sample_name_from_path(path):
    """Derive the sample name from a flagstat path.

    Mapping outputs are 05_mapping/{sample}/{sample}.flagstat, so the basename with the
    extension stripped is the sample name.
    """
    return os.path.basename(path).rsplit('.flagstat', 1)[0]


def main():
    parser = argparse.ArgumentParser(
        description="Compile samtools flagstat reports (default format) into a single "
                    "wide TSV, one row per sample."
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='One or more samtools flagstat files'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output TSV file path'
    )

    args = parser.parse_args()

    print(f"Found {len(args.input)} flagstat file(s) to process", file=sys.stderr)

    failed = []
    processed_count = 0

    with open(args.output, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(HEADER)

        for path in sorted(args.input):
            sample_name = sample_name_from_path(path)
            stats = parse_flagstat(path)

            if stats is not None:
                writer.writerow([sample_name] + [stats.get(col, 'NA') for col in HEADER[1:]])
                processed_count += 1
                print(f"Processed: {sample_name} ({stats['status']})", file=sys.stderr)
            else:
                failed.append((sample_name, path))

    print(f"\nCompleted: {processed_count}/{len(args.input)} samples processed",
          file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)

    # A non-empty but unparseable flagstat means the mapping step did not finish for that
    # sample. Fail loudly rather than shipping a summary that is silently missing rows.
    if failed:
        print(f"\nERROR: {len(failed)} flagstat report(s) could not be parsed:",
              file=sys.stderr)
        for sample_name, path in failed:
            print(f"  - {sample_name}: {path}", file=sys.stderr)
        print("These samples are absent from the summary. Check the corresponding "
              "reference_mapping logs; a truncated flagstat usually means bwa mem or "
              "samtools did not finish.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
