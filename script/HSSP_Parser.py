import pandas as pd
import numpy as np
import re
import pickle


def matcher(df, pattern):
    # The lines from which sections start
    lines = df[df.line.str.contains(pattern, regex=True)]
    return lines


def fetching_header(header_df):
    """
    :param header_df: header section
    :return: header dict keys: 'VERSION', 'PDBID', 'SEQLENGTH', 'NCHAIN', 'NALIGN'
    """
    header = {}
    # important line in header section
    # HSSP, PDBID, SEQLENGTH, NCHAIN, NALIGN
    pattern = 'HSSP|PDBID|SEQLENGTH|NCHAIN|NALIGN'
    lines = matcher(header_df, pattern).values.tolist()
    header['VERSION'] = lines[0][0].split('VERSION')[1].lstrip()
    header['PDBID'] = lines[1][0].split('PDBID')[1].strip()
    header['SEQLENGTH'] = int(lines[2][0].split('SEQLENGTH')[1].strip())
    header['NCHAIN'] = int(lines[3][0].split('NCHAIN')[1].split("chain")[0].strip())
    header['NALIGN'] = int(lines[4][0].split('NALIGN')[1].strip())

    return header


def remove_white_space(df):
    """
    :param df: dataframe
    :return: sec_list list
    """

    pattern = re.compile(r'\s+')
    sec_list = df.line.apply(lambda x: re.sub(pattern, ',', x.lstrip().rstrip()).split(',')).values.tolist()

    return sec_list


def fetching_proteins(proteins_sec):
    """
    :param proteins_sec: dataframe of protein section
    :return:proteins_dict dict
    """

    first_part = proteins_sec.line.apply(lambda x: x[:90])
    description_part = proteins_sec.line.apply(lambda x: x[90:])

    # remove whitespaces
    proteins_sec_list = remove_white_space(pd.DataFrame(first_part))
    for i, line_list in enumerate(proteins_sec_list):
        if i == 0:
            line_list[0] = line_list[0][:-1]
            line_list.append('PROTEIN')
        else:
            del line_list[1]
            if len(line_list) != 12:
                line_list.insert(2, ' ')
            line_list.append(description_part.iloc[i])

    protein_sec_df = pd.DataFrame(proteins_sec_list[1:], columns=['NR', 'ID', 'STRID', 'IDE', 'WSIM', 'IFIR', 'ILAS',
                                                                  'JFIR', 'JLAS', 'LALI', 'NGAP', 'LGAP', 'LSEQ2',
                                                                  'ACCNUM',  'PROTEIN'])
    # convert to int64
    protein_sec_df = protein_sec_df.astype({'NR': int, 'IDE': float, 'WSIM': float, 'IFIR': int, 'ILAS': int,
                                            'JFIR': int, 'JLAS': int, 'LALI': int, 'NGAP': int, 'LGAP': int,
                                            'LSEQ2': int})

    # set AliNo as index number
    protein_sec_df.set_index('NR', inplace=True)
    proteins_dict = {'PROTEINS': protein_sec_df.to_dict('index')}

    return proteins_dict


def fetching_alignment(alignment_sec, line_number):
    """
    :param alignment_sec: dataframe of alignment section
    :param line_number: list of line numbers which contain of '## ALIGNMENTS'
    :return:alignment_dict dict
    """

    # remove header lines
    first_line_this_part = line_number[0]
    removed_rows = []
    for ndx in line_number:
        removed_rows.extend([ndx - first_line_this_part, ndx - first_line_this_part + 1])
    alignment_sec = alignment_sec.drop(alignment_sec.index[removed_rows]).reset_index()

    analysis_part_alignment_section = alignment_sec.line.apply(lambda x: x[:51] + x[121:])
    # remove duplicates
    analysis_part_alignment_section.drop_duplicates(keep="first", inplace=True)
    # remove whitespaces
    analysis_part_alignment_sec_list = remove_white_space(pd.DataFrame(analysis_part_alignment_section))
    removed_ndx = []
    for ndx, line_list in enumerate(analysis_part_alignment_sec_list):
        line_list[1:3] = [''.join(line_list[1:3])]
        if line_list[-1].replace('.', '', 1).isnumeric():
            line_list.extend([np.nan, np.nan])
        if len(line_list) - 7 == 3:
            # no STRUCTURE
            line_list.insert(3, ' ')
        else:
            line_list[3:-7] = [''.join(line_list[3:-7])]
            line_list[3] = alignment_sec.line[ndx][17:26]
        if '!!' in line_list:
            removed_ndx.append(ndx)

    analysis_part_alignment_sec_list = [elem for i, elem in enumerate(analysis_part_alignment_sec_list)
                                        if i not in removed_ndx]
    # Putting alignments together
    alignment_part_alignment_section = alignment_sec.line.apply(lambda x: x[51:120])
    total_lines = len(analysis_part_alignment_sec_list)

    for i in range(total_lines):
        align = ''
        for j in range(len(line_number)):
            align = align + alignment_part_alignment_section[j * total_lines + i]
        analysis_part_alignment_sec_list[i].insert(-2, align)

    # create data frame
    align_sec_df = pd.DataFrame(analysis_part_alignment_sec_list[:], columns=['SeqNo', 'PDBNo', 'AA', 'STRUCTURE',
                                                                              'BP1', 'BP2', 'ACC', 'NOCC', 'VAR',
                                                                              'alignment', 'CHAIN', 'AUTHCHAIN'])

    # convert to int64
    align_sec_df = align_sec_df.astype({'SeqNo': int, 'ACC': int, 'NOCC': int, 'VAR': int})

    # set PDBNo as index number

    align_sec_df.set_index('PDBNo', inplace=True)
    alignment_dict = {'ALIGNMENTS': align_sec_df.to_dict('index')}

    return alignment_dict


def fetch_profile(profile_sec):
    """
    :param profile_sec: dataframe of profile section
    :return:profile_dict dict
    """
    profile_sec_list = remove_white_space(profile_sec)

    for line_list in profile_sec_list[1:]:
        line_list[1:3] = [''.join(line_list[1:3])]
        if line_list[-1].replace('.', '', 1).isnumeric():
            line_list.append(np.nan)
            line_list.append(np.nan)

    profile_sec_df = pd.DataFrame(profile_sec_list[1:], columns=profile_sec_list[0])

    # convert types
    for col in profile_sec_df.columns:
        if col not in ['PDBNo', 'CHAIN', 'AUTHCHAIN']:
            profile_sec_df[col] = pd.to_numeric(profile_sec_df[col])

    # set PDBNo as index number
    profile_sec_df.set_index('PDBNo', inplace=True)
    profile_dict = {'PROFILE': profile_sec_df.to_dict('index')}

    return profile_dict


def fetch_insertion_list(insertion_sec):
    """
    :param insertion_sec: dataframe of insertion section
    :return:insertion_dict dict keys: AliNo  IPOS  JPOS   Len Sequence
    """
    insertion_sec_list = remove_white_space(insertion_sec)
    removed_ndx = []
    for ndx, line_list in enumerate(insertion_sec_list):
        if '+' in line_list:
            insertion_sec_list[ndx - 1][-1] = insertion_sec_list[ndx - 1][-1] + line_list[-1]
            removed_ndx.append(ndx)
    insertion_sec_list = [elem for i, elem in enumerate(insertion_sec_list) if i not in removed_ndx]

    insertion_sec_df = pd.DataFrame(insertion_sec_list[1:-1], columns=insertion_sec_list[0])
    # convert to int64
    insertion_sec_df = insertion_sec_df.astype({'AliNo': int, 'IPOS': int, 'JPOS': int, 'Len': int})
    insertion_dict = {'INSERTION': insertion_sec_df.to_dict('index')}

    return insertion_dict


def parse_hssp(df):
    # The lines from which sections start
    pattern = "## PROTEINS : identifier and alignment statistics|## ALIGNMENTS|## SEQUENCE PROFILE AND ENTROPY|" \
              "## INSERTION LIST"
    lines_start_section = matcher(df, pattern)
    # The line number from which sections start
    lines_num_start_section = lines_start_section.index.values.tolist()

    # parsing sections
    hssp = {}
    # parsing header section
    hssp = fetching_header(df.loc[0:lines_num_start_section[0] - 1])
    # parsing PROTEINS : identifier and alignment statistics section
    hssp.update(fetching_proteins(df.loc[lines_num_start_section[0] + 1:lines_num_start_section[1] - 1].reset_index()))
    # parsing alignment section
    alignment_section_line_num = lines_num_start_section[1:-2]
    hssp.update(fetching_alignment(df.loc[lines_num_start_section[1]:lines_num_start_section[-2] - 1].reset_index(),
                                   alignment_section_line_num))
    # parsing SEQUENCE PROFILE AND ENTROPY
    hssp.update(fetch_profile(df.loc[lines_num_start_section[-2] + 1:lines_num_start_section[-1] - 1].reset_index()))
    # parsing INSERTION section
    hssp.update(fetch_insertion_list(df.loc[lines_num_start_section[-1] + 1:].reset_index()))

    return hssp


def write_pickle_file(hssp_dict, output_file):

    with open(output_file, 'wb') as handle:
        pickle.dump(hssp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_hssp_file(input_file, output_file):
    df = pd.read_table(input_file, names=['line'])

    # remove NOTATION
    df = df[~df.line.str.contains("NOTATION")]

    hssp_dict = parse_hssp(df)
    write_pickle_file(hssp_dict, output_file)

