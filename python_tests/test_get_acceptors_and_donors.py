from python_code.get_acceptors_and_donors import (
    dna_data_read,
    get_acceptors,
    get_donors,
    DONOR_SEQ,
    ACCEPTOR_SEQ,
)
import io
from python_tests.tests_utils import (
    TEST_DATA,
    SEQ_FROM_TEST_DATA,
    TEST_A,
    TEST_B,
    EXPECTED_TRUE_ACCEPTORS,
    EXPECTED_TRUE_DONORS,
    NO_OVERLAP_EXPECTED_FALSE_ACCEPTORS,
    NO_OVERLAP_EXPECTED_FALSE_DONORS,
    OVERLAP_EXPECTED_FALSE_ACCEPTORS,
    OVERLAP_EXPECTED_FALSE_DONORS
)


def test_dna_data_split():
    """
    Test of file is read in proper way
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    assert sequences == [SEQ_FROM_TEST_DATA]


def test_get_acceptors_dummy():
    """
    Dummy test which check if read sequences of length 2 ar only acceptor sequence.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_a, f_a = get_acceptors(0, 2, sequences, False)

    for x in t_a + f_a:
        assert x == ACCEPTOR_SEQ


def test_get_acceptors_dummy_2():
    """
    Check if acceptors sequence is in correct position.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    a = 2
    b = 3
    t_a, f_a = get_acceptors(a, b, sequences, False)

    for x in t_a + f_a:
        assert x[a : a + 2] == ACCEPTOR_SEQ


def test_get_donors_dummy():
    """
    Dummy test which check if read sequences of length 2 ar only donor sequence.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_d, f_d = get_donors(0, 2, sequences, False)

    for x in t_d + f_d:
        assert x == DONOR_SEQ


def test_get_donors_dummy_2():
    """
    Check if donor sequence is in correct position.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    a = 2
    b = 3
    t_a, f_a = get_donors(a, b, sequences, False)

    for x in t_a + f_a:
        assert x[a : a + 2] == DONOR_SEQ



def test_get_acceptors_no_overlap():
    """
    Test get acceptors for data mention in documentation, no overlap fragments.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_a, f_a = get_acceptors(TEST_A, TEST_B, sequences, False)

    for x in t_a:
        assert x in EXPECTED_TRUE_ACCEPTORS

    for x in f_a:
        assert x in NO_OVERLAP_EXPECTED_FALSE_ACCEPTORS


def test_get_acceptors_overlap():
    """
    Test get acceptors for data mention in documentation, overlap fragments.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_a, f_a = get_acceptors(TEST_A, TEST_B, sequences, True)

    for x in EXPECTED_TRUE_ACCEPTORS:
        assert x in t_a

    for x in OVERLAP_EXPECTED_FALSE_ACCEPTORS:
        assert x in f_a


def test_get_donors_no_overlap():
    """
    Test get donors for data mention in documentation, no overlap fragments.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_d, f_d = get_donors(TEST_A, TEST_B, sequences, False)

    for x in t_d:
        assert x in EXPECTED_TRUE_DONORS

    for x in f_d:
        assert x in NO_OVERLAP_EXPECTED_FALSE_DONORS


def test_get_donors_overlap():
    """
    Test get donors for data mention in documentation, overlap fragments.
    """
    with io.StringIO(TEST_DATA) as f:
        sequences = dna_data_read(f)

    t_d, f_d = get_donors(TEST_A, TEST_B, sequences, True)

    for x in t_d:
        assert x in EXPECTED_TRUE_DONORS

    for x in f_d:
        assert x in OVERLAP_EXPECTED_FALSE_DONORS
