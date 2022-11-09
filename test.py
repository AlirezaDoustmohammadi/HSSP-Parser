from script.HSSP_Parser import read_hssp_file
import pickle


def read_pickle_file(file_name):
    with open(file_name, 'rb') as handle:
        hssp_dict = pickle.load(handle)

    # 'VERSION', 'PDBID', 'SEQLENGTH', 'NCHAIN', 'NALIGN', 'PROTEINS', 'ALIGNMENTS', 'PROFILE', 'INSERTION'
    # explanation: 'VERSION': hssp version,
    # 'NCHAIN': number of chains
    # 'NALIGN': number of aligned proteins
    # 'PROTEINS' (dict): ## PROTEINS : identifier and alignment statistics --> key: NR
    # 'ALIGNMENTS' (dict): ## ALIGNMENTS --> key: PDBNo
    # 'PROFILE' (dict): ## SEQUENCE PROFILE AND ENTROPY --> key: PDBNo
    # 'INSERTION' (dict): ## INSERTION LIST --> key: index (start from zero)
    print(hssp_dict.keys())
    print(hssp_dict['PROTEINS'][1])
    print(hssp_dict['ALIGNMENTS']['438A'])


if __name__ == '__main__':
    if __name__ == '__main__':
        hssp_file = 'input file/1taq.hssp'
        output_file = 'output/1taq.pickle'
        read_hssp_file(hssp_file, output_file)
        # read output file
        read_pickle_file(output_file)
