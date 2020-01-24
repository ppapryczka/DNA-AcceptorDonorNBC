import argparse
import os
import re
import sys
from typing import List, Tuple

# sequence corresponding to donor
DONOR_SEQ: str = "GT"
# sequence corresponding to acceptor
ACCEPTOR_SEQ: str = "AG"

# default donors result file
DEFAULT_DONORS_RESULT_FILENAME: str = "donor.dat"
# default acceptors result file
DEFAULT_ACCEPTORS_RESULT_FILENAME: str = "acceptor.dat"
# default input data file
DEFAULT_DATA_FILE: str = "araclean.dat"

# default left value of length
DEFAULT_A_LEN: int = 10
# default right value of length
DEFAULT_B_LEN: int = 10


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
    parser.add_argument("A", help="'left' length of donor/acceptor", type=int)

    # B length argument
    parser.add_argument("B", help="'right' length of donor/acceptor", type=int)

    # input file argument
    parser.add_argument(
        "input",
        help="input file with sequences",
        type=lambda x: parser_check_if_file_exists(parser, x),
    )

    # donors result file argument
    parser.add_argument(
        "donors",
        help="donors result file",
        nargs="?",
        default=DEFAULT_DONORS_RESULT_FILENAME,
        type=str,
    )

    # acceptors results file argument
    parser.add_argument(
        "acceptors",
        help="acceptors result file",
        nargs="?",
        default=DEFAULT_ACCEPTORS_RESULT_FILENAME,
        type=str,
    )

    # parse arguments
    args = parser.parse_args(command_args)

    # read sequences
    with open(args.input) as dna_file:
        dna_sequences = dna_data_read(dna_file)

    # run get acceptors and donors
    true_donor, false_donor, true_acceptor, false_acceptor = get_acceptors_and_donors(args.A, args.B, dna_sequences)

    save_acceptor_and_donors_to_file(true_donor, false_donor, true_acceptor, false_acceptor, args.donors, args.acceptors)


def dna_data_read(file) -> List:
    """
    Open ``file``, read lines and take only
    ``Introns``, ``Exons`` and ``Data`` section of sequences.

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


def get_acceptors_and_donors(a_len: int, b_len: int, dna_sequences: List) -> Tuple[List, List, List, List]:
    """
    Read data from ``input_file``, get false and real donors and save them to``donors_output``, get false
    and real acceptors and save them to ``acceptors_output``.

    Args:
        a_len: Left length of donor/acceptor.
        b_len: Right length of donor/acceptor.
        dna_sequences: Sequences read from data file, each element of list is map which
        contains "Introns", "Exons", "Data".


    """
    # get donors
    true_donors: List = []
    false_donors: List = []
    for dna_sequence in dna_sequences:
        introns_begin_list = [x[0] for x in dna_sequence["Introns"]]
        true_donors.extend(
            get_true_fragments(introns_begin_list, a_len, b_len, dna_sequence["Sequence"])
        )
        false_donors.extend(
            get_false_fragments(
                introns_begin_list, a_len, b_len, dna_sequence["Sequence"], DONOR_SEQ
            )
        )

    # get acceptors
    true_acceptors: List = []
    false_acceptors: List = []
    for dna_sequence in dna_sequences:
        introns_end_list = [x[1] for x in dna_sequence["Introns"]]
        true_acceptors.extend(
            get_true_fragments(introns_end_list, a_len, b_len, dna_sequence["Sequence"])
        )
        false_acceptors.extend(
            get_false_fragments(
                introns_end_list, a_len, b_len, dna_sequence["Sequence"], ACCEPTOR_SEQ
            )
        )

    return true_donors, false_donors, true_acceptors, false_acceptors


def save_acceptor_and_donors_to_file(
    true_donors: List,
    false_donors: List,
    true_acceptors: List,
    false_acceptors: List,
    donors_output: str,
    acceptors_output: str,
):
    with open(donors_output, "w") as donors_f:
        for x in true_donors:
            donors_f.write(f"{1}\n{x}\n")
        for x in false_donors:
            donors_f.write(f"{0}\n{x}\n")

    with open(acceptors_output, "w") as acceptors_f:
        for x in true_acceptors:
            acceptors_f.write(f"{1}\n{x}\n")
        for x in false_acceptors:
            acceptors_f.write(f"{0}\n{x}\n")


def get_true_fragments(
    true_positions_list: List[int], a_len: int, b_len: int, sequence_data: str
):
    """
    Get real donors/acceptors from given ``sequence_data`` described by ``true_fragments_positions``.
    """
    true_fragments = []
    for x in true_positions_list:
        if not is_fragment_outside_of_sequence(a_len, b_len, x, len(sequence_data)):
            true_fragments.append(sequence_data[x - a_len: x + b_len])
    return true_fragments


def get_false_fragments(
    true_fragments_positions: List[int], a_len: int, b_len: int, sequence_data: str, fragment: str
):
    """
    Get false donors/acceptors from given ``sequence_data`` which have ``fragment at position ``A``. 
    ``true_fragments_positions`` are use to omit real donors/acceptors.
    """
    false_fragments = []
    possible_false_fragments = [m.start() for m in re.finditer(fragment, sequence_data)]
    for false_fragment_pos in possible_false_fragments:
        if is_good_false_fragment(
            true_fragments_positions, a_len, b_len, false_fragment_pos, len(sequence_data)
        ):
            false_fragments.append(
                sequence_data[false_fragment_pos - a_len: false_fragment_pos + b_len]
            )
    return false_fragments


def line_of_numbers_to_tuples(line: str) -> List[Tuple[int, int]]:
    """
    Convert line of number split be white signs to list of tuples.
    
    Returns:
        List of tuples. 
    """
    result = []
    number_splited = line.split()

    number_iter = iter(number_splited)
    for x in number_iter:
        result.append((int(x), int(next(number_iter))))

    return result


def is_fragment_outside_of_sequence(a_len: int, b_len: int, pos: int, seq_len: int) -> bool:
    """
    Check if pos give a valid fragment which not go outside sequence.
    """
    if pos + a_len >= seq_len:
        return True
    elif pos - b_len < 0:
        return True

    return False


def is_good_false_fragment(
    true_positions_list: List, a_len: int, b_len: int, pos: int, seq_len: int
) -> bool:
    """
    Check if fragment given by ``pos`` is not outside sequence and don't overlap any of true donors/acceptors
    mentioned in ``true_positions_list``.
    """
    if is_fragment_outside_of_sequence(a_len, b_len, pos, seq_len):
        return False

    if have_collision_with_true_fragment(true_positions_list, a_len, b_len, pos):
        return False

    return True


def have_collision_with_true_fragment(
    true_positions_list: List, a_len: int, b_len: int, pos: int
) -> bool:
    """
    Check if false fragment given by it ``position`` don't have overlapping with given 
    true fragments in ``introns_list``.
    """
    for true_pos in true_positions_list:
        if true_pos - a_len <= pos - a_len <= true_pos + b_len:
            return True

        if true_pos - a_len <= pos + b_len <= true_pos + b_len:
            return True

    return False


if __name__ == "__main__":
    get_acceptors_and_donors_command(sys.argv[1:])
