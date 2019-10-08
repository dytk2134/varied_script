#!/usr/bin/env python
# Contributed by Li-Mei Chiang <dytk2134 [at] gmail [dot] com> (2019)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.


import sys
import logging
import subprocess
import os
import yaml
from multiprocessing import Pool
from tqdm import tqdm

__version__ = '1.0.0'

# logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)


def config_parser(config_file):
    if not os.path.exists(config_file):
        logger.error('Failed to find the config file: %s' % (config_file))
        sys.exit(1)

    with open(config_file, 'r') as c_in:
        try:
            config_dict = yaml.safe_load(c_in)
        except yaml.YAMLError as exc:
            logger.error(exc)
            sys.exit(1)

    return config_dict

def allele_freq(chromosome_set, variant_dict, variant_list, db, db_path, populations):
    if not os.path.exists(db_path):
        logger.error('%s: No such file or directory' % (db_path))
        sys.exit(1)

    not_found = set(range(len(variant_list)))
    for chromosome in chromosome_set:
        ref_data = os.path.join(db_path, chromosome)
        if os.path.exists(ref_data):
            with open(ref_data, 'r') as in_f:
                for line in in_f:
                    line = line.strip()
                    if line:
                        # [chr, pos, ref, alt, ref_freq, alt_freq]
                        tokens = line.split('\t')
                        key = (tokens[0], tokens[1], tokens[2], tokens[3])
                        if key in variant_dict:
                            for idx in variant_dict[key]:
                                if idx not in not_found:
                                    continue
                                else:
                                    not_found.remove(idx)
                                    variant_list[idx].extend(tokens[4:])
    # not found
    for idx in not_found:
        tmp = ['.'] * (len(populations) * 2)
        variant_list[idx].extend(tmp)
    return variant_list
def annotation(chromosome_set, variant_dict, variant_list, anno_table, db_path):
    if not os.path.exists(db_path):
        logger.error('%s: No such file or directory' % (db_path))
        sys.exit(1)
    not_found = set(range(len(variant_list)))
    for chromosome in chromosome_set:
        ref_data = os.path.join(db_path, chromosome)
        if os.path.exists(ref_data):
            with open(ref_data, 'r') as in_f:
                for line in in_f:
                    line = line.strip()
                    if line:
                        # [chr, pos, ref, alt, id]
                        tokens = line.split('\t')
                        key = (tokens[0], tokens[1], tokens[2], tokens[3])
                        if key in variant_dict:
                            for idx in variant_dict[key]:
                                if idx not in not_found:
                                    continue
                                else:
                                    not_found.remove(idx)
                                    variant_list[idx].extend(tokens[4])
    # not found
    for idx in not_found:
        variant_list[idx].append('.')
    return variant_list

def annovar(variant_dict, variant_list, header, config_dict, output_dir, annovar_input, annovar_output):
    try:
        if config_dict['Tools']['annovar']['include'] == True:
            opperation = list()
            protocol = list()
            if config_dict['Tools']['annovar']['gene_based_annotation']:
                opperation.extend(['g'] * len(config_dict['Tools']['annovar']['gene_based_annotation']))
                protocol.extend(config_dict['Tools']['annovar']['gene_based_annotation'])
            if config_dict['Tools']['annovar']['region_based_annotation']:
                opperation.extend(['r'] * len(config_dict['Tools']['annovar']['region_based_annotation']))
                protocol.extend(config_dict['Tools']['annovar']['region_based_annotation'])
            if config_dict['Tools']['annovar']['filter_based_annotation']:
                opperation.extend(['f'] * len(config_dict['Tools']['annovar']['filter_based_annotation']))
                protocol.extend(config_dict['Tools']['annovar']['filter_based_annotation'])

            cmd = [
                os.path.join(config_dict['Tools']['annovar']['tool_path'], 'table_annovar.pl'),
                annovar_input,
                config_dict['Tools']['annovar']['humandb_path'],
                '-buildver',
                config_dict['Tools']['annovar']['buildver'],
                '-remove',
                '--thread',
                str(config_dict['Common']['threads']),
                '-protocol',
                ','.join(protocol),
                '-operation',
                ','.join(opperation),
                '-nastring',
                '.',
                '--outfile',
                annovar_output
            ]

            subprocess.call(' '.join(cmd), shell=True)

            annovar_result = annovar_output + '.hg19_multianno.txt'
            annovar_header = list()
            annovar_items_count = 0
            not_found = set(range(len(variant_list)))
            if os.path.exists(annovar_result):
                with open(annovar_result, 'r') as in_f:
                    for line in in_f:
                        line = line.strip()
                        if line:
                            tokens = line.split('\t')
                            if not annovar_header:
                                annovar_header = tokens
                                header.extend(tokens[5:])
                                annovar_items_count = len(tokens[5:])
                            else:
                                key = (tokens[0], tokens[1], tokens[3], tokens[4])
                                if key in variant_dict:
                                    for idx in variant_dict[key]:
                                        if idx not in not_found:
                                            continue
                                        else:
                                            not_found.remove(idx)
                                            variant_list[idx].extend(tokens[5:])
                for idx in not_found:
                    tmp = ['.'] * (annovar_items_count)
                    variant_list[idx].extend(tmp)
    except KeyError:
        pass
    return variant_list, header
def main(input_files, config_file, output_dir):
    # read config file (yaml format)
    config_dict = config_parser(config_file)

    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except IOError as err:
            logger.error(err)
            sys.exit(1)

    for input_file in tqdm(list(input_files), desc='Annotating ... '):
        variant_dict = dict()
        variant_list = list()
        chromosome_set = set()
        header = ['chromosome', 'position', 'alt', 'ref']
        annovar_filename = os.path.splitext(os.path.basename(input_file))[0] + '.avinput'
        with open(os.path.join(output_dir, annovar_filename), 'w') as out_f:
            with open(input_file, 'r') as in_f:
                count = 0
                for line in in_f:
                    line = line.strip()
                    if line:
                        tokens = line.split('\t')
                        chromosome = tokens[0].upper().replace('CHR', '')
                        chromosome_set.add(chromosome)
                        key = (chromosome, tokens[1], tokens[2], tokens[3])
                        if key not in variant_dict:
                            variant_dict[key] = list()

                        variant_dict[key].append(count)
                        variant_list.append([chromosome, tokens[1], tokens[2], tokens[3]])
                        avinput_output = [chromosome, tokens[1], str(int(tokens[1]) + len(tokens[2]) - 1), tokens[2], tokens[3]]
                        out_f.write('\t'.join(avinput_output)+'\n')
                    count += 1
            # annovar
            annovar_input = os.path.join(output_dir, annovar_filename)
            annovar_output = os.path.join(output_dir, os.path.splitext(os.path.basename(input_file))[0])
            variant_list, header = annovar(variant_dict, variant_list, header, config_dict, output_dir, annovar_input, annovar_output)

            # database
            try:
                # allele freq
                for freq_table in config_dict['Databases']['allele_freq']:
                    if config_dict['Databases']['allele_freq'][freq_table]['include'] == True:
                        for pop in config_dict['Databases']['allele_freq'][freq_table]['populations']:
                            header.extend(['%s_%s_Ref' % (freq_table, pop), '%s_%s_Alt' % (freq_table, pop)])
                        variant_list = allele_freq(chromosome_set, variant_dict, variant_list, freq_table, config_dict['Databases']['allele_freq'][freq_table]['db_path'], config_dict['Databases']['allele_freq'][freq_table]['populations'])
            except KeyError:
                pass
            # annotation
            try:
                for anno_table in config_dict['Databases']['annotation']:
                    if config_dict['Databases']['annotation'][anno_table]['include'] == True:
                        header.append(anno_table)
                        variant_list = annotation(chromosome_set, variant_dict, variant_list, anno_table, config_dict['Databases']['annotation'][anno_table]['db_path'])
            except KeyError:
                pass

        with open(os.path.join(output_dir, os.path.splitext(os.path.basename(input_file))[0] + '_varied.tsv'), 'w') as out_f:
            out_f.write('\t'.join(header) + '\n')
            for varinat in variant_list:
                out_f.write('\t'.join(varinat) + '\n')

if __name__ == '__main__':
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\

    Quick start:
    %(prog)s -in_b input.bigwig -regions region.txt -f genome.fasta -out_b output.bigwig
    """))
    # argument
    parser.add_argument('-i', '--input_files', nargs='+', help='Input tsv file')
    parser.add_argument('-c', '--conf', help='conf file (yaml format)')
    parser.add_argument('-d', '--output_dir', type=str, help='Specify the directory name for saving results. default: varied_output', default='varied_output')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit(1)

    if not args.input_files:
        logger.error('Please specify an input file with the argument "-i".')
        sys.exit(1)

    if not args.conf:
        logger.error('Please specifiy the config file with the argument "-c".')
        sys.exit(1)

    main(input_files=args.input_files, config_file=args.conf, output_dir=args.output_dir)