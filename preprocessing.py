import os
import sys
import re
import argparse
import subprocess
import logging
from misc import check_file_exists, check_create_dir, execute_subprocess

logger = logging.getLogger()


"""
=============================================================
HEADER
=============================================================

INSTITUTION:IiSGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
VERSION=0.1
CREATED: 21 Sep 2020
REVISION: 

================================================================
END_OF_HEADER
================================================================
"""


def fastqc_quality(r1, r2, output_dir, threads=8):
    check_create_dir(output_dir)

    cmd = ['fastqc', r1, r2, '-o', output_dir, '--threads', str(threads)]

    execute_subprocess(cmd)


def fastp_trimming(r1, r2, sample, output_dir, threads=6, min_qual=20, window_size=10, min_len=35):
    check_create_dir(output_dir)

    output_trimmed_r1 = os.path.join(
        output_dir, sample + ".trimmed_R1.fastq.gz")
    output_trimmed_r2 = os.path.join(
        output_dir, sample + ".trimmed_R2.fastq.gz")

    html_dir = os.path.join(output_dir, 'html')
    json_dir = os.path.join(output_dir, 'json')

    check_create_dir(html_dir)
    check_create_dir(json_dir)

    html_file = os.path.join(html_dir, sample + '_fastp.html')
    json_file = os.path.join(json_dir, sample + '_fastp.json')

    cmd = ['fastp',
           '--in1', r1,
           '--in2', r2,
           '--out1', output_trimmed_r1,
           '--out2', output_trimmed_r2,
           '--detect_adapter_for_pe',
           '--cut_tail',
           '--cut_window_size', str(window_size),
           '--cut_mean_quality', str(min_qual),
           '--length_required', str(min_len),
           '--json', json_file,
           '--html', html_file,
           '--thread', str(threads)]

    #    '--trim_front1', '3',
    #    '--trim_front2', '3',
    #    '--trim_tail1', '3',
    #    '--trim_tail2', '3',

    execute_subprocess(cmd)

# FORMAT TO COMPARE FASTQC OUTPUT


def extract_sample_from_fastq(R_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R = os.path.basename(R_file)

    long_suffix = re.search('_S.*', basename_R)
    dot_suffix = re.search('.R[12].*', basename_R)

    if long_suffix:
        match = long_suffix.group()
        sample_name = basename_R.split(match)[0]
    elif dot_suffix:
        match = dot_suffix.group()
        sample_name = basename_R.split(match)[0]
    else:
        sample_name = basename_R

    return str(sample_name)


def extract_processed_html(fastqc_folder, sample, read_type):
    for root, _, files in os.walk(fastqc_folder):
        for name in files:
            fileName = os.path.join(root, name)
            if ('Quality/processed' in fileName) and (read_type in name) and (name.endswith('fastqc.html') and (name.startswith(sample))):
                return fileName


def extract_files_html(fastqc_folder):
    html_pairs = {}
    count = 0
    for root, _, files in os.walk(fastqc_folder):
        for name in files:
            fileName = os.path.join(root, name)
            if 'Quality/raw' in fileName and name.endswith('fastqc.html'):
                sample = extract_sample_from_fastq(fileName)
                if "R1" in name:
                    read_type = 'R1'
                elif "R2" in name:
                    read_type = 'R2'
                raw_html = fileName
                processed_html = extract_processed_html(
                    fastqc_folder, sample, read_type)
                html_pairs[count] = [raw_html, processed_html]
                count = count + 1
    return html_pairs


def extract_quality_graph(html_file):
    with open(html_file, 'r') as f:
        content = f.read()
        image_tag = re.search(
            r'<img class="indented" src=.*alt="Per base quality graph" .+?/>', content)
    return image_tag.group(0)


def extract_basic_stats(html_file):
    with open(html_file, 'r') as f:
        content = f.read()
        table_tag = re.search(
            r'Basic Statistics</h2>(<table>(.+?)</table>)', content)

    return table_tag.group(1)


def format_html_image(output_folder):
    files = extract_files_html(output_folder)
    html_template = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>covidma quality output</title>
    <style type="text/css">
        body {
        margin: 0 auto;
        }
        </style>
    </head>
    <body>
    IMAGEHTMLPAIRED
    </body>
    </html>
    """
    output_file = os.path.join(output_folder, 'fastq_image_report.html')
    all_images_tables = ''
    for number, pair in files.items():
        logger.debug(number)
        logger.debug(pair)
        div_structure = """
        <div class="container">
        <table>
            <tr>
            <th>FILENAMERAW</th>
            <th>FILENAMETRIMMED</th>
            </tr>
            <tr>
            <td>TABLEQUALRAW</td>
            <td>TABLEQUALTRIMMED</td>
            </tr>
            <tr>
            <td>IMAGEQUALRAW</td>
            <td>IMAGEQUALTRIMMED</td>
            </tr>
        </table>
        <br>
        </div>
        """
        div_structure = div_structure.replace('FILENAMERAW', pair[0])
        div_structure = div_structure.replace('FILENAMETRIMMED', pair[1])
        div_structure = div_structure.replace(
            'TABLEQUALRAW', extract_basic_stats(pair[0]))
        div_structure = div_structure.replace(
            'TABLEQUALTRIMMED', extract_basic_stats(pair[1]))
        div_structure = div_structure.replace(
            'IMAGEQUALRAW', extract_quality_graph(pair[0]))
        div_structure = div_structure.replace(
            'IMAGEQUALTRIMMED', extract_quality_graph(pair[1]))

        all_images_tables = all_images_tables + div_structure

    final_html_template = html_template.replace(
        'IMAGEHTMLPAIRED', all_images_tables)
    with open(output_file, 'w+') as f:
        f.write(final_html_template)
