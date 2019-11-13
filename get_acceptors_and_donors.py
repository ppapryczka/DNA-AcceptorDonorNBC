from typing import List,  Tuple,  Dict
import re
import argparse
import sys
import os

DONOR_SEQ: str = "GT"
ACCEPTOR_SEQ: str = "AG"

DEFAULT_DONORS_RESULT_FILENAME = "donor.dat"
DEFAULT_ACCEPTORS_RESULT_FILENAME = "acceptor.dat"


def get_acceptors_and_donors_command(command_args: List[str]) -> None:
    """
    Create parser, parse args given in ``command_args`` and run ``get_acceptors_and_donors``. 
    Command take four arguments:  A, B, input data, result donors output file, 
    result acceptors output file.

    Args:
        command_args: Arguments for command.
    """    
    
    # create parser 
    parser = argparse.ArgumentParser()
    
    # A length argument
    parser.add_argument(
        "A",  
        help ="'left' length of donor/acceptor", 
        type=int
    )
    
    # B length argument
    parser.add_argument(
        "B", 
        help="'right' length of donor/acceptor", 
        type=int
    )
    
    # input file argument
    parser.add_argument(
        "input", 
        help="input file with sequences",  
        type=lambda x: parser_check_if_file_exists(parser,  x)
    )
    
    # donors result file argument 
    parser.add_argument(
        "donors",  
        help="donors result file",  
        nargs = "?", 
        default= DEFAULT_DONORS_RESULT_FILENAME,  
        type=str
    )
    
    # acceptors results file argument
    parser.add_argument(
        "acceptors",  
        help="acceptors result file",  
        nargs = "?", 
        default= DEFAULT_ACCEPTORS_RESULT_FILENAME,  
        type=str
    )
    
    # parse arguments
    args = parser.parse_args(command_args)
    
    # run get acceptors and donors
    get_acceptors_and_donors(args.A,  args.B,  args.input,  args.donors,  args.acceptors)
    
    
def parser_check_if_file_exists(parser: argparse.ArgumentParser, file_path: str) -> str:
    """
    Check if file exists - if not call parser.error, else return ``file_path``
    as a result.

    Args:
        parser: Command parser.
        file_path: Path to file.
    """
    if os.path.isfile(file_path):
        return file_path
    else:
        parser.error(f"The file {file_path} doesn't exist!")
    
def get_acceptors_and_donors(A: int,  B: int,  input_file: str,  donors_output: str,  acceptors_output: str) -> None:
    """
    Read data from ``input_file``, get false and real donors and save them to``donors_output``, get false
    and real acceptors and save them to ``acceptors_output``.
    """
    
    dna_sequences = dna_data_split(input_file)
    
    # get donors 
    true_d= []
    false_d= []
    for dna_sequence in dna_sequences:
        introns_begin_list = [ x[0] for x in dna_sequence["Introns"] ]
        true_d.extend(get_true_fragments(introns_begin_list,  A,  B,  dna_sequence["Sequence"]))
        false_d.extend(get_false_fragments(introns_begin_list,  A,  B,  dna_sequence["Sequence"],  DONOR_SEQ))

    # get acceptors 
    true_a= []
    false_a= []
    for dna_sequence in dna_sequences:
        introns_end_list = [ x[1] for x in dna_sequence["Introns"] ]
        true_a.extend(get_true_fragments(introns_end_list,  A,  B,  dna_sequence["Sequence"]))
        false_a.extend(get_false_fragments(introns_end_list,  A,  B,  dna_sequence["Sequence"],  ACCEPTOR_SEQ))
    
    with open(donors_output,  "w") as donors_f:
        for x in true_d:
            donors_f.write(f"{1}\n{x}\n")
        for x in false_d:
            donors_f.write(f"{0}\n{x}\n")
    
    with open(acceptors_output,  "w") as acceptors_f:
        for x in true_a:
            acceptors_f.write(f"{1}\n{x}\n")
        for x in false_a:
            acceptors_f.write(f"{0}\n{x}\n")
    
    print("We are on a mission from God!")
    
def dna_data_split(data_file: str ) -> List:
    """
    Open ``data_file``, read lines and take only 
    ``Introns``, ``Exons`` and ``Data`` section of files. 
    
    Returns:
        List of dictionaries, each have following structure: 
        {"Introns", "Exons", "Sequence"}. 
    """
    # read lines
    with open(data_file) as dna_f:
        dna_file_lines: List[str] =  [line.rstrip('\n') for line in dna_f]
    
    lines_iter = iter(dna_file_lines)
    
    # read data
    dna_sequences = []
    dna_sequence = {}
    for line in lines_iter:
        # introns section
        if line == "Introns":
            dna_sequence["Introns"] = line_of_numbers_to_tuples(next(lines_iter))
        # exons section
        elif line == "Exons":
            dna_sequence["Exons"] = line_of_numbers_to_tuples(next(lines_iter))
        # data section
        elif line == "Data":
            dna_sequence["Sequence"] = next(lines_iter)
            # save sequence
            dna_sequences.append(dna_sequence)
            dna_sequence = {}
    
    return dna_sequences

def get_true_fragments(true_positions_list: List[int],  A: int,  B: int, sequence_data: str):
    """
    Get real donors/acceptors from given ``sequence_data`` described by ``true_fragments_positions``.
    """
    true_fragments = []
    for x in true_positions_list:
        if not is_fragment_outside_of_sequence(A,  B,  x,  len(sequence_data)):
            true_fragments.append(sequence_data[x-A:x+B])
    return true_fragments

def get_false_fragments(true_fragments_positions: int,  A: int,  B: int, sequence_data: str, fragment: str):
    """
    Get false donors/acceptors from given ``sequence_data`` which have ``fragment at position ``A``. 
    ``true_fragments_positions`` are use to omit real donors/acceptors.
    """
    false_fragments = []
    possible_false_fragments = [m.start() for m in re.finditer(fragment,  sequence_data)]
    for false_fragment_pos in possible_false_fragments:
            if  is_good_false_fragment(true_fragments_positions,  A,  B,  false_fragment_pos,  len(sequence_data)):
                false_fragments.append(sequence_data[false_fragment_pos-A:false_fragment_pos+B])
    return false_fragments

def line_of_numbers_to_tuples(line: str) -> List[Tuple[int,  int]]:
    """
    Convert line of number split be white signs to list of tuples.
    
    Returns:
        List of tuples. 
    """
    result = []
    number_splited=  line.split()
    
    number_iter = iter(number_splited)
    for x in number_iter:
        result.append((int(x),  int(next(number_iter))))
    
    return result

def is_fragment_outside_of_sequence(A: int,  B: int,  pos: int,  seq_len: int) -> bool:
    """
    Check if pos give a valid fragment which not go outside  sequence. 
    """
    if pos + B >= seq_len:
       return True
    elif pos - A < 0:
       return True
       
    return False

def is_good_false_fragment(true_positions_list: int,  A:int,  B: int, pos: int,  seq_len: int) -> bool:
    """
    Check if fragment given by ``pos`` is not outside sequence and don't overlap any of true donors/acceptors
    mentioned in ``true_positions_list``.
    """
    if  is_fragment_outside_of_sequence(A,  B,  pos,  seq_len):
        return False
    
    if  have_coolision_with_true_fragment(true_positions_list,  A,  B,  pos):
        return False
    
    return True

def have_coolision_with_true_fragment(true_positions_list: int,  A:int,  B: int,  pos: int ) ->bool:
    """
    Check if false fragment given by it ``position`` don't have overlapping with given 
    true fragments in ``introns_list``.
    """
    for true_pos in true_positions_list:
        if true_pos-A <= pos-A <= true_pos+B:
            return True
        
        if true_pos-A <= pos+B <= true_pos+B:
            return True

    return False

if __name__ == "__main__":
    get_acceptors_and_donors_command(sys.argv[1:])
