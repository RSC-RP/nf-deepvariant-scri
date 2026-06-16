#! /usr/bin/env python3

import csv

# Script to take input CSV of BAMs and convert it to a pedigree for Somalier

def read_lookup(file):
    '''
    Read a text file with sample IDs from the metadata sheet in the first column,
    and sample IDs from the BAM header in the second column. Convert to dict.
    '''
    out = dict()
    with open(file, mode = 'rt', newline = '') as incon:
        rdr = csv.reader(incon, delimiter = '\t')
        for row in rdr:
            out[row[0]] = row[1]
    return out

def read_samplesheet(file):
    '''
    Read the sample sheet that was input to the pipeline, and extract information
    needed for building pedigree.
    '''
    out = list()
    keepkeys = ['family', 'proband_sex', 'proband_id', 'father_id', 'mother_id']
    with open(file, mode = 'rt', newline = '') as incon:
        rdr = csv.DictReader(incon)
        for row in rdr:
            out.append({k: row[k] for k in keepkeys})
    return out

def assemble_pedigree(fams, lookup):
    '''
    Generate a list of pedigree rows to write.
    '''
    # Make preliminary pedigree rows from the file.
    out = set()
    for row in fams:
        dad = row['father_id']
        if dad == '':
            dad = '-9'
        else:
            dad = lookup[dad]
        mom = row['mother_id']
        if mom == '':
            mom = '-9'
        else:
            mom = lookup[mom]
        if row['proband_sex'] in {'Male', 'male', 'M'}:
            sex = '1'
        else:
            assert row['proband_sex'] in {'Female', 'female', 'F'}
            sex = '2'

        proband_row = (
            row['family'],
            lookup[row['proband_id']],
            dad,
            mom,
            sex,
            '-9'
        )
        out.add(proband_row)

        if dad != '-9' and dad not in {row[1] for row in out}:
            out.add((
                row['family'],
                dad,
                '-9',
                '-9',
                '1',
                '-9'
            ))
        if mom != '-9' and mom not in {row[1] for row in out}:
            out.add((
                row['family'],
                mom,
                '-9',
                '-9',
                '2',
                '-9'
            ))
    # check for any duplicate individuals, keep those with parents
    all_ind = {row[1] for row in out}
    out2 = set()
    for ind in all_ind:
        theserow = {row for row in out if row[1] == ind}
        if len(theserow) > 1:
            fams = {row[0] for row in theserow}
            dads = {row[2] for row in theserow if row[2] != '-9'}
            moms = {row[3] for row in theserow if row[3] != '-9'}
            sex = {row[4] for row in theserow}
            if(len(dads) > 1 or len(moms) > 1 or len(sex) > 1):
                raise Exception("Problem building pedigree.")
            if len(dads) == 1:
                dad = dads.pop()
            else:
                dad = '-9'
            if len(moms) == 1:
                mom = moms.pop()
            else:
                mom = '-9'
            out2.add((fams.pop(), ind, dad, mom, sex.pop(), '-9'))
        else:
            out2.update(theserow)
    # consolidate family names
    out3 = set()
    while len(out2) > 0:
        theserows = {out2.pop()}
        search = True
        # Gather up all individuals in this family
        while search:
            famind = {row[1] for row in theserows} | \
                {row[2] for row in theserows if row[2] != '-9'} | \
                {row[3] for row in theserows if row[3] != '-9'}
            newrows = {row for row in out2 if row[1] in famind or row[2] in famind or row[3] in famind}
            if len(newrows) == 0:
                search = False
            else:
                theserows.update(newrows)
                out2.difference_update(newrows)
        # Find all family names and pick one
        fams = {row[0] for row in theserows}
        fam = fams.pop()
        outrows = {(fam, row[1], row[2], row[3], row[4], row[5]) for row in theserows}
        out3.update(outrows)
    return out3

# Testing family consolidation and renaming
#out = {
#    ('fam1', 'C', 'B', 'A', '1', '-9'),
#    ('fam1', 'A', '-9', '-9', '2', '-9'),
#    ('fam1', 'B', '-9', '-9', '1', '-9'),
#    ('fam2', 'B', 'D', 'E', '1', '-9'),
#    ('fam2', 'D', '-9', '-9', '1', '-9'),
#    ('fam2', 'E', '-9', '-9', '2', '-9')
#    }

if __name__ == "__main__":
    import sys
    bam_csv = sys.argv[1]
    sample_lookup = sys.argv[2]
    out_ped = sys.argv[3]
    sample_dict = read_lookup(sample_lookup)
    families = read_samplesheet(bam_csv)
    ped_rows = assemble_pedigree(families, sample_dict)
    with open(out_ped, mode = 'wt', newline = '') as outcon:
        wtr = csv.writer(outcon, delimiter = '\t', quoting = csv.QUOTE_NONE)
        for row in ped_rows:
            wtr.writerow(row)
