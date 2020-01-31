import argparse
import os
import re
import sys
from typing import List, Tuple
from enum import Enum

# sequence corresponding to donor
DONOR_SEQ: str = "GT"
# sequence corresponding to acceptor
ACCEPTOR_SEQ: str = "AG"

# default donors result file
DEFAULT_RESULT_OUTPUT_FILENAME: str = "result.dat"
# default input data file
DEFAULT_INPUT_DATA_FILE: str = "araclean.dat"

# default left value of length
DEFAULT_A_LEN: int = 10
# default right value of length
DEFAULT_B_LEN: int = 10


class DnaFragmentType(Enum):
    DONOR = 1
    ACCEPTOR = 2


def parser_check_if_file_exists(parser: argparse.ArgumentParser, file_path: str) -> str:
    """
    Check if file exists - if not call parser.error, else return ``file_path`` as a result.

    Args:
        parser: Command parser.
        file_path: Path to file.

    Returns:
        Original ``file_path`` when file exists.
    """
    if os.path.isfile(file_path):
        return file_path
    else:
        parser.error(f"The file {file_path} doesn't exist!")


def parser_check_if_given_type_is_correct(parser: argparse.ArgumentParser, seq_type: str) -> str:
    """
    Check if type is supported - if not call parser.error, else return ``seq_type`` as a result.

    Args:
        parser: Command parser.
        seq_type: Type of sequence - look ``DnaFragmentType``.

    Returns:
        Original ``seq_type`` when file exists.
    """
    if seq_type in [f'{dft.name}' for dft in DnaFragmentType]:
        return seq_type
    else:
        parser.error(f"Value {seq_type} don't supported! Supported types: {[f'{dft.name}' for dft in DnaFragmentType]}")


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
    parser.add_argument("-A", "--a_len", help="'left' length of donor/acceptor", required=True, type=int)

    # B length argument
    parser.add_argument("-B", "--b_len", help="'right' length of donor/acceptor", required=True, type=int)

    # overlap
    parser.add_argument(
        "-o",
        "--overlap",
        action="store_true",
        default = False,
        help="if flag use generate addition false donors/acceptors with overlap fragments with true donors/acceptors",
    )

    # input file argument
    parser.add_argument(
        "-i",
        "--input",
        help="input file with sequences",
        default=DEFAULT_INPUT_DATA_FILE,
        required=True,
        type=lambda x: parser_check_if_file_exists(parser, x),
    )

    # type
    parser.add_argument(
        "-t",
        "--type",
        help=f"type of DNA fragment Supported types: {[f'{dft.name}' for dft in DnaFragmentType]}",
        required = True,
        type= lambda x: parser_check_if_given_type_is_correct(parser, x)
    )

    # result
    parser.add_argument(
        "-r",
        "--result",
        help="result file",
        default=DEFAULT_RESULT_OUTPUT_FILENAME,
        type=str,
        required=True
    )

    # parse arguments
    args = parser.parse_args(command_args)

    # read sequences
    with open(args.input) as dna_file:
        dna_sequences = dna_data_read(dna_file)

    # generate output
    if args.type == DnaFragmentType.ACCEPTOR.name:
        true_acceptor, false_acceptor = get_acceptors(args.a_len, args.b_len, dna_sequences, args.overlap)
        save_sequences_to_file(true_acceptor, false_acceptor, args.result)
    elif args.type == DnaFragmentType.DONOR.name:
        true_donors, false_donors = get_donors(args.a_len, args.b_len, dna_sequences, args.overlap)
        save_sequences_to_file(true_donors, false_donors, args.result)
    else:
        raise Exception("Wrong type given!")


def dna_data_read(file) -> List:
    """
    Open ``file``, read lines and take only ``Introns``, ``Exons`` and ``Data`` section of sequences.

    Args:
        file: Open file like object with sequences.

    Returns:
        List of dictionaries, each have following structure:
        {"Introns", "Exons", "Sequence"}.
    """
    # read lines
    dna_file_lines: List[str] = [line.rstrip("\n") for line in file]

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


def get_acceptors(
    a_len: int, b_len: int, dna_sequences: List, overlap_fragments: bool
) -> Tuple[List, List]:
    """
    Get false and real acceptors from ``dna_sequences``.

    Args:
        a_len: Left length of acceptor.
        b_len: Right length of acceptor.
        dna_sequences: Sequences read from data file, each element of list is map which
        contains "Introns", "Exons", "Data".
        overlap_fragments: If true generate overlap fragments, if false don't generate overlap fragments.

    Returns:
        True acceptors, false acceptors lists.
    """
    # get acceptors
    true_acceptors: List = []
    false_acceptors: List = []
    for dna_sequence in dna_sequences:

        # IMPORTANT -1!!!
        introns_end_list = [x[1] - 1 for x in dna_sequence["Introns"]]
        true_acceptors.extend(
            get_true_fragments(introns_end_list, a_len, b_len, dna_sequence["Sequence"])
        )
        false_acceptors.extend(
            get_false_fragments(
                introns_end_list, a_len, b_len, dna_sequence["Sequence"], ACCEPTOR_SEQ, overlap_fragments
            )
        )

    return true_acceptors, false_acceptors


def get_donors(
    a_len: int, b_len: int, dna_sequences: List, overlap_fragments: bool
) -> Tuple[List, List]:
    """
    Get false and real donors from ``dna_sequences``.

    Args:
        a_len: Left length of donor.
        b_len: Right length of donor.
        dna_sequences: Sequences read from data file, each element of list is map which
        contains "Introns", "Exons", "Data".
        overlap_fragments: If true generate overlap fragments, if false don't generate overlap fragments.

    Returns:
        True donors, false donors lists.
    """
    # get donors
    true_donors: List = []
    false_donors: List = []
    for id, dna_sequence in enumerate(dna_sequences):
        introns_begin_list = [x[0] for x in dna_sequence["Introns"]]
        true_donors.extend(
            get_true_fragments(
                introns_begin_list, a_len, b_len, dna_sequence["Sequence"]
            )
        )
        false_donors.extend(
            get_false_fragments(
                introns_begin_list, a_len, b_len, dna_sequence["Sequence"], DONOR_SEQ, overlap_fragments
            )
        )

    return true_donors, false_donors


def save_sequences_to_file(
    true_seq: List,
    false_seq: List,
    output: str
) -> None:
    """
    Save acceptors/donors to given files.

    Args:
        true_seq: List of true donors/acceptors.
        false_seq: List of false donors/acceptors.
        output: Output file name.
    """
    with open(output, "w") as seq_f:
        for x in true_seq:
            if x != "":
                seq_f.write(f"{1}\n{x}\n")
        for x in false_seq:
            if x != "":
                seq_f.write(f"{0}\n{x}\n")

def get_true_fragments(
    true_positions_list: List[int], a_len: int, b_len: int, sequence_data: str
):
    """
    Get real donors/acceptors from given ``sequence_data`` described by ``true_fragments_positions``.

    Args:
        true_positions_list: List of true positions of donors/acceptors.
        a_len: Left length of donor/acceptor.
        b_len: Right length of donor/acceptor.
        sequence_data: DNA data.

    Returns:
        List of true fragments.
    """
    true_fragments = []
    for x in true_positions_list:
        if not is_fragment_outside_of_sequence(a_len, b_len, x, len(sequence_data)):
            true_fragments.append(sequence_data[x - a_len : x + b_len])
    return true_fragments


def get_false_fragments(
    true_fragments_positions: List[int],
    a_len: int,
    b_len: int,
    sequence_data: str,
    fragment: str,
    overlap_fragments: bool
):
    """
    Get false donors/acceptors from given ``sequence_data`` which have ``fragment at position ``A``. 
    ``true_fragments_positions`` are use to omit real donors/acceptors.

    Args:
        true_fragments_positions: List of true positions of donors/acceptors.
        a_len: Left length of donor/acceptor.
        b_len: Right length of donor/acceptor.
        sequence_data: DNA data.
        fragment: Fragment searched in ``sequence_data`` - look ACCEPTOR_SEQ, DONOR_SEQ.
        overlap_fragments: If true generate overlap fragments, if false don't generate overlap fragments.

    Returns:
        List of false fragments.
    """
    false_fragments = []
    possible_false_fragments = [m.start() for m in re.finditer(fragment, sequence_data)]

    for false_fragment_pos in possible_false_fragments:
        if is_good_false_fragment(
            true_fragments_positions,
            a_len,
            b_len,
            false_fragment_pos,
            len(sequence_data),
            overlap_fragments
        ):
            false_fragments.append(
                sequence_data[false_fragment_pos - a_len : false_fragment_pos + b_len]
            )
    return false_fragments


def line_of_numbers_to_tuples(line: str) -> List[Tuple[int, int]]:
    """
    Convert line of number split be white signs to list of tuples.

    Args:
        line: Line to split numbers from.

    Returns:
        List of tuples. 
    """
    result = []
    number_split = line.split()

    number_iter = iter(number_split)
    for x in number_iter:
        result.append((int(x), int(next(number_iter))))

    return result


def is_fragment_outside_of_sequence(
    a_len: int, b_len: int, pos: int, seq_len: int
) -> bool:
    """
    Check if pos give a valid fragment which not go outside sequence.

    Args:
        a_len: Left length of donor/acceptor.
        b_len: Right length od donor/acceptor.
        pos: position to check.
        seq_len: Length of original DNA sequence.

    Returns:
        True if fragment goes outside sequence, false if is correct.
    """
    if pos + a_len >= seq_len:
        return True
    elif pos - b_len < 0:
        return True

    return False


def is_good_false_fragment(
    true_positions_list: List, a_len: int, b_len: int, pos: int, seq_len: int, overlap_fragments: bool
) -> bool:
    """
    Check if fragment given by ``pos`` is not outside sequence and don't overlap any of true donors/acceptors
    mentioned in ``true_positions_list``.

    Args:
        true_positions_list: List of true positions of donors/acceptors.
        a_len: Left length of donor/acceptor.
        b_len: Right length od donor/acceptor.
        pos: Position to check.
        seq_len: Length of original DNA sequence.
        overlap_fragments: If true accept overlap fragments, if false don't accept overlap fragments.

    Returns:
        True if fragment is 'good'(don't overlap true fragment and don't goes outside sequence), false instead.
    """
    if is_fragment_outside_of_sequence(a_len, b_len, pos, seq_len):
        return False

    if have_collision_with_true_fragment(true_positions_list, a_len, b_len, pos, overlap_fragments):
        return False

    return True


def have_collision_with_true_fragment(
    true_positions_list: List, a_len: int, b_len: int, pos: int, overlap_fragments: bool
) -> bool:
    """
    Check if false fragment given by it ``position`` don't have overlapping with given 
    true fragments in ``true_positions_list``.

    Args:
        true_positions_list: List of true positions of donors/acceptors.
        a_len: Left length of donor/acceptor.
        b_len: Right length od donor/acceptor.
        pos: Position to check.
        overlap_fragments: If true accept overlap fragments, if false don't accept overlap fragments.

    Returns:
        True if fragment have collision with true fragment, false if don't have collision.
    """
    if not overlap_fragments:
        for true_pos in true_positions_list:
            if true_pos - a_len <= pos - a_len <= true_pos + b_len:
                return True

            if true_pos - a_len <= pos + b_len <= true_pos + b_len:
                return True
    else:
        for true_pos in true_positions_list:
            if true_pos == pos:
                return True

    return False


if __name__ == "__main__":
    get_acceptors_and_donors_command(sys.argv[1:])
