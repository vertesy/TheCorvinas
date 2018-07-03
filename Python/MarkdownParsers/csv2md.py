#!/usr/bin/env python
# Source: https://gist.github.com/ghing/ba983b90b137872fa7b5

import csv
import sys

import click
"You need to install click !!!!"


def underlines(fieldnames):
    underlines = []
    for fieldname in fieldnames:
        underline = ""
        for i in range(len(fieldname)):
            underline += "-"

        underlines.append(underline)

    return underlines

def order_row(fieldnames, row):
    return [row[fn] for fn in fieldnames]


@click.command()
@click.argument('input', type=click.File('rb'), default=None, required=False)
def csv2md(input):
    """
    Output a table formatted with GitHub Flavored Markdown from CSV
    
    Examples:

    Output Markdown for the first 8 rows of an Excel file.

         openelex/us/ia/cache/20100608__ia__primary__adair__precinct.xls | head -n 8 | ./csv2md.py
    
    """
    if input is None:
        input = sys.stdin

    reader = csv.DictReader(input)
    sys.stdout.write(' | '.join(reader.fieldnames) + "\n")
    sys.stdout.write(' | '.join(underlines(reader.fieldnames)) + "\n")
    for row in reader:
        sys.stdout.write(' | '.join(order_row(reader.fieldnames, row)) + "\n")

if __name__ == '__main__':
    print 11
    csv2md()