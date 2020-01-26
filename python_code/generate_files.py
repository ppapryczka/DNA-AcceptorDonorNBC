from python_code.classification_NBC import get_data_to_classifier
from python_code.get_acceptors_and_donors import (
    DEFAULT_DATA_FILE,
    dna_data_split,
    get_acceptors_and_donors,
    DEFAULT_A,
    DEFAULT_B,
)

if __name__ == "__main__":
    with open(DEFAULT_DATA_FILE) as f:
        sequences = dna_data_split(f)

    a = 40
    b = 40

    true_d, false_d, true_a, false_a = get_acceptors_and_donors(a, b, sequences)
    df = get_data_to_classifier(true_d, false_d, a, b)

    a = 15
    b = 15

    true_d, false_d, true_a, false_a = get_acceptors_and_donors(a, b, sequences)
    df = get_data_to_classifier(true_d, false_d, a, b)

    a = 35
    b = 35

    true_d, false_d, true_a, false_a = get_acceptors_and_donors(a, b, sequences)
    df = get_data_to_classifier(true_d, false_d, a, b)

    """
    for a in range(5, 21):
        for b in range(5, 21):
            true_d, false_d, true_a, false_a = get_acceptors_and_donors(a, b, sequences)
            df = get_data_to_classifier(true_d, false_d, a, b)
    """
