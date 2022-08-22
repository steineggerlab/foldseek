#include "MSANode.h"
#include "Matcher.h"
#include "Parameters.h"

#include <cassert>

const char* binary_name = "test_generatetree";


// Gap in query: no Q match state for on T match state (special case for speed reasons)
// i:       i-1   i-1       di=0
// Q:  XXXXXX.....~~~XXX
// T:  YYYYYYyyyyyYYYYYY
// j:       j-1   j
// l:       lprev l         dl=6
void test_dq_0() {
    std::string newSeq  = "AAAA";
    std::string targetSeq = "AAAABBBBCCCC";
    std::cout << "Testing dq == 0\n";
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;

    newSeq.append(noQueryMS(targetSeq, 3, 8));
    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "AAAAbbbbc");
    std::cout << std::endl;
}

// One Q match state for one T match state (treated as special case for speed reasons)
// i:       i-1   i         di=1
// Q:  XXXXXX.....XXXXXX
// T:  YYYYYYyyyyyYYYYYY
// j:       j-1   j
// l:       lprev l         dl=6
void test_dq_1() {
    std::string newSeq = "AAAA";
    std::string targetSeq = "AAAABBBBCCCC";
    std::cout << "Testing dq == 1\n";
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;

    newSeq.append(oneQueryMS(targetSeq, 3, 8));
    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "AAAAbbbbC");
    std::cout << std::endl;
}

// More Match states in Q than Inserts in the T sequence
// => half T inserts y left, half right-aligned in uc, gaps to fill up
// Number of T insert residues to be left-aligned: (int)(dl/2)
// i:        iprev  i       di=7
// Q:  XXXXXXXXXXXXXXXXXX
// T:  YYYYYYYyyy-yyYYYYY
// j:        j-1    j
// l:        lprev  l       dl=6
void test_dq_gt_dt() {
    //                Query: XXXXXXXXXXXXXXXXXX
    std::string querySeq =  "XXXXXXXXXXXXXXXXXX";
    std::string newSeq    = "XXXXXXX";
    std::string targetSeq = "XXXXXXXyyyyyXXXXX";
    std::cout << "Testing dq >= dt\n";
    std::cout << "   querySeq = " << querySeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    gtQueryMS(newSeq, targetSeq, 6, 12, 7, 6);

    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "XXXXXXXYYY-YYX");
    std::cout << std::endl;

    newSeq    = "XXXXXXX";
    targetSeq = "XXXXXXXyy-yyXXXXX";
    std::cout << "Testing with gap\n";
    std::cout << "   querySeq = " << querySeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    gtQueryMS(newSeq, targetSeq, 6, 12, 7, 6);
    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "XXXXXXXYY--YYX");

    std::cout << std::endl;
}

// Fewer Match states in Q than inserts in T sequence
// => half of available space di for left
// - half for right-aligned T inserts,
// rest in lc
// number of T inserts to be left-aligned in uc: (int)(di/2),
// i:        iprev i       di=5
// Q:  XXXXXXXXX.XXXXXXX
// T:  YYYYYYYyyyyyYYYYY
// j:        j-1   j
// l:        lprev l       dl=6
void test_dq_lt_dt() {
    std::string newSeq    = "XXXXXXX";
    std::string targetSeq = "XXXXXXXyyyyyXXXXX";
    std::cout << "Testing dq < dt\n";
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;
    ltQueryMS(newSeq, targetSeq, 6, 12, 5, 6);
    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "XXXXXXXYYyYYX");

    std::cout << "\nTesting with gap\n";
    newSeq    = "XXXXXXX";
    targetSeq = "XXXXXXXyy-yyXXXXX";
    std::cout << " newSeq (b) = " << newSeq << std::endl;
    std::cout << "  targetSeq = " << targetSeq << std::endl;
    ltQueryMS(newSeq, targetSeq, 6, 12, 5, 6);
    std::cout << " newSeq (a) = " << newSeq << std::endl;
    assert(newSeq == "XXXXXXXYYYYX");
    std::cout << std::endl;   
}

void test_full_MSA() {
    std::vector<MSASequence> seqs;
    
    MSASequence seq0(0, "ABCDEFGHIJKLMNOPQRSTUVWXYZ\n");
    MSASequence seq1(1, "EFFGHIJKLMNOPQRSTUVW\n");
    MSASequence seq2(2, "CDEGHIJKLMNOPQRSTUV\n");
    MSASequence seq3(3, "GHIJKLMMNOPQRSTUVWXYZ\n");
    MSASequence seq4(4, "GHIJKLQRSTUVWXYZ\n");
    MSASequence seq5(5, "GHIJKLMMNOPMNOPQRSTUVWXYZ\n");
    MSASequence seq6(6, "ABCDEFFFLMNOPQRSTUVWXYZ\n");

    seq0.print();
    seq1.print();
    seq2.print();
    seq3.print();
    seq4.print();
    seq5.print();
    seq6.print();

    // After iterating best hits and creating MSANodes
    std::vector<MSASequence *> node_one_seqs = { &seq0 };
    std::vector<MSASequence *> node_two_seqs = { &seq1 };
    std::vector<MSASequence *> node_thr_seqs = { &seq2 };
    std::vector<MSASequence *> node_fou_seqs = { &seq3 };
    std::vector<MSASequence *> node_fiv_seqs = { &seq4 };
    std::vector<MSASequence *> node_six_seqs = { &seq5 };
    std::vector<MSASequence *> node_sev_seqs = { &seq6 };

    MSANode node_one(0, node_one_seqs);
    MSANode node_two(1, node_two_seqs);
    MSANode node_thr(2, node_thr_seqs);
    MSANode node_fou(3, node_fou_seqs);
    MSANode node_fiv(4, node_fiv_seqs);
    MSANode node_six(5, node_six_seqs);
    MSANode node_sev(6, node_sev_seqs);

    // Test cases
    Matcher::result_t result = {};

    // Query ID=0, Target ID=1
    //       
    //        4                  22
    // 0: ABCDEF.GHIJKLMNOPQRSTUVWXYZ
    //        MMDMMMMMMMMMMMMMMMMM
    // 1:     EFFGHIJKLMNOPQRSTUVW
    //        0                  19
    //
    result.dbKey = 1;
    result.backtrace = "MMDMMMMMMMMMMMMMMMMM";
    result.qStartPos = 4;
    result.qEndPos = 22;
    result.dbStartPos = 0;
    result.dbEndPos = 19;
    node_one.addNode(0, &node_two, result);


    // Fewer match states in Q than inserts in T
    // Q:  XXXXXXXXX.XXXXXXX
    // T:  YYYYYYYyyyyyYYYYY
    //
    // 0: ABCDEFGHIJKMNOPQRSTUVWXYZ
    // 6: ABCDEFGHIJKLLLLLMNOPQRSTUVWXYZ
    //
    /* alignNewSequence( */
    /*     "XXXXXXXXX.XXXXXXX", */
    /*     "YYYYYYYyyyyyYYYYY", */
    /*     0, 0, 16, */
    /*     result, */
    /*     matches */
    /* ) */


    // Query ID=0, Target ID=2
    //
    //      2                 21
    // 0: ABCDEFGHIJKLMNOPQRSTUVWXYZ
    //
    //        MMDMMMMMMMMMMMMMMMMM
    // 1:     EFFGHIJKLMNOPQRSTUVW
    //
    //      MMMIMMMMMMMMMMMMMMMM
    // 2:   CDE-GHIJKLMNOPQRSTUV
    //      0                  18
    //
    result.dbKey = 2;
    result.backtrace = "MMMIMMMMMMMMMMMMMMMM";
    result.qStartPos = 2;
    result.qEndPos = 21;
    result.dbStartPos = 0;
    result.dbEndPos = 18;
    node_one.addNode(0,  &node_thr, result);

    // Single insertions
    // Query ID=3, Target ID=2
    //
    //       0                   20
    // 3:    GHIJKLMMNOPQRSTUVWXYZ
    //       MMMMMMMIMMMMMMMMM
    // 2: CDEGHIJKLM-NOPQRSTUV
    //       3               18
    //
    result.dbKey = 2;
    result.backtrace = "MMMMMMMIMMMMMMMMM";
    result.qStartPos = 0;
    result.qEndPos = 20;
    result.dbStartPos = 3;
    result.dbEndPos = 18;
    node_fou.addNode(3, &node_one, result);
    /* return EXIT_SUCCESS; */


    // Multiple insertions
    // Query ID=3, Target ID=4
    //
    //    0                   21
    // 3: GHIJKLMMNOPQRSTUVWXYZ
    //    MMMMMMIIIIIMMMMMMMMMM
    // 4: GHIJKL-----QRSTUVWXYZ
    //    0                   16
    //
    result.dbKey = 4;
    result.backtrace = "MMMMMMIIIIIMMMMMMMMMM";
    result.qStartPos = 0;
    result.qEndPos = 21;
    result.dbStartPos = 0;
    result.dbEndPos = 16;
    node_fou.addNode(3, &node_fiv, result);
    /* return EXIT_SUCCESS; */

    // Multiple deletions
    // Query ID=3, Target ID=5
    //
    //    0                       21
    // 3: GHIJKLMMNOP....QRSTUVWXYZ
    //    MMMMMMMMMMMDDDDMMMMMMMMMM
    // 5: GHIJKLMMNOPMNOPQRSTUVWXYZ
    //    0                       25
    //
    result.dbKey = 5;
    result.backtrace = "MMMMMMMMMMMDDDDMMMMMMMMMM";
    result.qStartPos = 0;
    result.qEndPos = 21;
    result.dbStartPos = 0;
    result.dbEndPos = 25;
    node_fou.addNode(3, &node_six, result);

    // Multiple insertions and deletions
    // Query ID=3, Target ID=5
    // 
    //    0                          26
    // 0: ABCDEF--GHIJKLMNOPQRSTUVWXYZ
    //    MMMMMMDDIIIIIMMMMMMMMMMMMMMM
    // 6: ABCDEFFF-----LMNOPQRSTUVWXYZ
    //    0                          23
    //
    result.dbKey = 6;
    result.backtrace = "MMMMMMDDIIIIIMMMMMMMMMMMMMMM";
    result.qStartPos = 0;
    result.qEndPos = 26;
    result.dbStartPos = 0;
    result.dbEndPos = 23;
    node_fou.addNode(0, &node_sev, result);
}

std::string mapAndAlign(std::string query, std::string target, Matcher::result_t result) {
    size_t qId = findAlnStartIdx(query, result.qStartPos);
    size_t tId = findAlnStartIdx(target, result.dbStartPos);
    std::cout << "\nqId: " << qId << ", tId: " << tId << std::endl;
    std::vector<int> matches = mapBacktrace(
        query,
        target,
        qId,
        tId,
        result.backtrace
    );

    std::cout << "Matches: \n";
    for (int member : matches) std::cout << member << " ";
    std::cout << std::endl;

    std::string newSeq = alignNewSequence(
        query,
        target,
        qId,
        tId,
        query.length(),
        result,
        matches
    );
    return newSeq;
}

void testAndPrint(std::string query, std::string target, std::string newSeq, std::string cigar, std::string expected) {
    std::cout << "     Query: " << query << std::endl;
    std::cout << "     CIGAR: " << cigar << std::endl;
    std::cout << "    Target: " << target << std::endl;
    std::cout << "  Expected: " << expected << std::endl;
    std::cout << "    newSeq: " << newSeq << std::endl;
    assert(newSeq == expected + '\n');
}

void test_alignNewSequence() {
    std::cout << "Testing alignNewSequence()\n";

    Matcher::result_t result;
    std::string querySeq;
    std::string targetSeq; 
    std::string newSeq;
    std::string expected;

    // All matches
    querySeq         = "ABCDEFGHIJKL";
    targetSeq        = "ABCDEFGHIJKL";
    expected         = "ABCDEFGHIJKL";
    result.backtrace = "MMMMMMMMMMMM";
    result.qStartPos = 0;
    result.qEndPos = 11;
    result.dbStartPos = 0;
    result.dbEndPos = 11;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    std::cout << "Testing all match:\n";
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);

    // Insertion
    querySeq         = "ABCDEFGHIJKL";
    targetSeq        = "ABCDEGHIJKL";
    expected         = "ABCDE-GHIJKL";
    result.backtrace = "MMMMMIMMMMMM";
    result.dbEndPos = 10;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    std::cout << "Testing insertion:\n";
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);

    // Deletion
    querySeq         = "ABCDEGHIJKL";
    targetSeq        = "ABCDEFGHIJKL";
    expected         = "ABCDEfGHIJKL";
    result.backtrace = "MMMMMDMMMMMM";
    result.qEndPos = 10;
    result.dbEndPos = 11;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);

    // Insertion and deletion
    querySeq         = "ABCDEGHIJKL";
    targetSeq        = "ABCDEFGHJKL";
    expected         = "ABCDEfGH-JKL";
    result.backtrace = "MMMMMDMMIMMM";
    result.qEndPos = 10;
    result.dbEndPos = 10;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    std::cout << "Testing both:\n";
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);

    // Starts with an insertion
    querySeq         = "ABCDEFGHIJKL";
    targetSeq        = "aBCDEFGHIJKL";
    expected         = "aBCDEFGHIJKL";
    result.backtrace = "MMMMMMMMMMMM";
    result.qEndPos = 11;
    result.dbEndPos = 11;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    std::cout << "Testing insertion start:\n";
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);

    // Already aligned
    querySeq         = "BCDEFGhiJKL";
    targetSeq        = "ABCdeFGHIJKL";
    expected         = "ABCDEFGHIJKL";
    result.backtrace = "DMMMMMMMMMMM";
    result.qEndPos = 10;
    result.dbEndPos = 11;
    newSeq = mapAndAlign(querySeq, targetSeq, result);
    std::cout << "Testing pre-aligned:\n";
    testAndPrint(querySeq, targetSeq, newSeq, result.backtrace, expected);
}

int main(int, const char**) {
    test_dq_0();
    test_dq_1();
    test_dq_gt_dt();
    test_dq_lt_dt();
    test_alignNewSequence();
    return EXIT_SUCCESS;

    test_full_MSA();
    return EXIT_SUCCESS;
}
