from typing import Dict

TEST_DATA: str = \
"""
>Seq 21 Len:
37
Introns
 4 27
Exons
 0 3 28 36
Data
TTTTGTCAGCACATACGTAGCGAGTTCAGATGTGCT
"""

SEQ_FROM_TEST_DATA: Dict = \
{
    "Introns":[(4, 27)],
    "Exons": [(0, 3), (28, 36)],
    "Sequence" : "TTTTGTCAGCACATACGTAGCGAGTTCAGATGTGCT"
}