import unittest
from bam_slice.bam_slice import merge_overlap_intervals, handle_nones, get_subsequence_indexes

class TestMergeOverlapIntervals(unittest.TestCase):

    def test_NoOverlappingIntervals_NoChange(self):
        intervals = [
            [2, 4],
            [6, 9],
            [11, 12],
        ]

        result = merge_overlap_intervals(intervals)
        expected = [
            [2, 4],
            [6, 9],
            [11, 12],
        ]
        self.assertEqual(result, expected)

    def test_TwoIntervalsEqualEndStart_NoChange(self):
        intervals = [
                [6, 9],
                [9, 12]
                ]

        result = merge_overlap_intervals(intervals)
        expected = [
                [6, 9],
                [9, 12]
                ]
        self.assertEqual(result, expected)

    def test_TwoIntervalsOverlap_Merge(self):
        intervals = [
                [6, 9],
                [8, 12]
                ]

        result = merge_overlap_intervals(intervals)
        expected = [
                [6, 12]
                ]
        self.assertEqual(result, expected)

    def test_ThreeIntervalsOverlap_Merge(self):
        intervals = [
                [6, 9],
                [8, 12],
                [11, 14]
                ]

        result = merge_overlap_intervals(intervals)
        expected = [
                [6, 14]
                ]
        self.assertEqual(result, expected)

    def test_ThreeIntervalsOverlapTwoEqualsEndStart_MergeOverlapDontMergeEquals(self):
        intervals = [
                [6, 9],
                [8, 12],
                [11, 14],
                [14, 16]
                ]

        result = merge_overlap_intervals(intervals)
        expected = [
                [6, 14],
                [14, 16]
                ]
        self.assertEqual(result, expected)


class TestHandleNones(unittest.TestCase):
    def test_NoNones_NoChange(self):
        aligned_pairs = [
                (1, 1),
                (2, 2)
                ]

        result = handle_nones(aligned_pairs)
        expected = ([1, 2], [1, 2])
        self.assertEqual(result, expected)

    def test_OneNoneRef_ChangeToPreviousRefPosition(self):
        aligned_pairs = [
                (1, 1),
                (2, None)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([1, 2], [1, 1])
        self.assertEqual(result, expected)

    def test_NoneAtRefStart_ChangeStartToMinusOne(self):
        aligned_pairs = [
                (1, None),
                (2, 1)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([1, 2], [-1, 1])
        self.assertEqual(result, expected)

    def test_AllRefNones_AllRefMinusOnes(self):
        aligned_pairs = [
                (1, None),
                (2, None)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([1, 2], [-1, -1])
        self.assertEqual(result, expected)

    def test_OneNoneRead_ChangeToPreviousReadPosition(self):
        aligned_pairs = [
                (1, 1),
                (None, 2)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([1, 1], [1, 2])
        self.assertEqual(result, expected)

    def test_NoneAtReadStart_ChangeStartToMinusOne(self):
        aligned_pairs = [
                (None, 0),
                (0, 1)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([-1, 0], [0, 1])
        self.assertEqual(result, expected)

    def test_AllReadNones_AllReadMinusOnes(self):
        aligned_pairs = [
                (None, 1),
                (None, 2)
                ]
        result = handle_nones(aligned_pairs)
        expected = ([-1, -1], [1, 2])
        self.assertEqual(result, expected)


class TestGetSubsequenceIndexes(unittest.TestCase):
    def test_NoRepeatedPositions_SimpleOneToOneMapping(self):
        interval = [2, 5]
        read_positions = list(range(30, 36))
        ref_positions = list(range(1, 7))
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (31, 34)
        self.assertEqual(result, expected)

    def test_IntervalAtStartOfRef_ZeroInStartIndex(self):
        interval = [0, 2]
        read_positions = list(range(30, 36))
        ref_positions = [-1, -1, -1, 0, 0, 1]
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (33, 35)
        self.assertEqual(result, expected)

    def test_RepeatedRefPositionsSomeInInterval_GetLeftAndRightOfRepeatsInInterval(self):
        interval = [2, 5]
        read_positions = list(range(30, 36))
        ref_positions = [1, 1, 1, 2, 2, 2]
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (33, 35)
        self.assertEqual(result, expected)

    def test_AllRefPositionsBeforeStartInterval_None(self):
        interval = [2, 5]
        read_positions = list(range(30, 36))
        ref_positions = [-1, -1, -1, 0, 0, 1]
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = None
        self.assertEqual(result, expected)

    def test_RefPositionsRepeatOnEndInterval_AllReadPositions(self):
        interval = [2, 5]
        read_positions = list(range(30, 36))
        ref_positions = [4, 5, 5, 5, 5, 5]
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (30, 35)
        self.assertEqual(result, expected)

    def test_ReadPositionsRepeatAcrossInterval_ReadPostionTwice(self):
        interval = [2, 5]
        read_positions = [29, 30, 30, 30, 30, 31]
        ref_positions = list(range(1, 7))
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (30, 30)
        self.assertEqual(result, expected)

    def test_TwoReadIndexesRepeating_BothReadIndexes(self):
        interval = [2, 5]
        read_positions = [30, 30, 30, 31, 31, 31]
        ref_positions = list(range(1, 7))
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (30, 31)
        self.assertEqual(result, expected)

    def test_ReadIndexInIntervalContainsMinusOnes_ChangeMinusOnesToZero(self):
        interval = [2, 5]
        read_positions = [-1, -1, -1, 0, 0, 1]
        ref_positions = list(range(1, 7))
        clean_aligned_pairs = (read_positions, ref_positions)
        result = get_subsequence_indexes(interval, clean_aligned_pairs)
        expected = (0, 0)
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
