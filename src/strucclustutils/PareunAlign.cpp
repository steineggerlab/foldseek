//
// Created by Martin Steinegger on 1/15/21.
//
//
//  Created by Stephanie Kim on 4/28/21.
//

#include "PareunAlign.h"
#include "EvalueComputation.h"

#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;


// creating an array that stores both integers and string.
struct alignmentInfo {
    int score;
    unsigned int query_start;
    unsigned int query_end;
    unsigned int target_start;
    unsigned int target_end;
    string cigar_string;
};

// Score matrix for adding gap alignments
vector<vector<int>> initializeScoreMatrix(Sequence & qSeq, Sequence & tSeq, SubstitutionMatrix &subMat) {

    vector<vector<int>> alignMatrix(qSeq.L, vector<int> (tSeq.L,0));
    vector<int> Scores;

    for (int i = 0; i < (qSeq.L)-1; i++) {
        for (int j = 0; j < (tSeq.L)-1; j++) {

            // Find the max score by computing 4 options:
            Scores.push_back(subMat.subMatrix[qSeq.numSequence[i]][tSeq.numSequence[j]]);
            Scores.push_back(subMat.subMatrix[qSeq.numSequence[i]][tSeq.numSequence[j+1]]);
            Scores.push_back(subMat.subMatrix[qSeq.numSequence[i+1]][tSeq.numSequence[j]]);
            Scores.push_back(0);

            // Put max score to alignMatrix (filling the align matrix)
            sort(Scores.begin(), Scores.end());
            alignMatrix[i+1][j+1] = Scores.back();
            Scores.clear();
        }
    }
    return alignMatrix;
}

// finding top 4 ungapped local alignments from alignMatrix. (score,row,col,length)
vector<vector<int>> ungappedAlign (Sequence & qSeq, Sequence & tSeq, vector<vector<int>> & A) {

    vector<vector<int>> maxScorelist {{}}; // return vector
    // Going through columns (Target sequence)
    for (int m = 0; m < (qSeq.L)-1; m++) {

        int MaxScore = 0; int Score = 0; int m2 = m; int n = 0;
        int alnLen=0; int TotalnLen=0; int endrow=0; int endcol=0;

        while (m2 < (qSeq.L)-1 && n < (tSeq.L)-1) {

            if (tSeq.getSeqData()[n] == qSeq.getSeqData()[m2]) {
                Score += A[m2+1][n+1];
                alnLen++;
            }

            else {
                if (Score > MaxScore) {
                    MaxScore = Score; Score = 0;
                    TotalnLen = alnLen; alnLen = 0;
                    endrow = n+1; endcol = m2+1;
                }
            }
            m2++; n++;
        }

        MaxScore=Score; TotalnLen=alnLen; endrow=m2+1; endcol=n+1; // necessary to get the last values

        if (MaxScore != 0) {
            maxScorelist.push_back({MaxScore, TotalnLen, endrow, endcol});
        }
    }
    // Going through rows (Query sequence)
    for (int m = 0; m < (tSeq.L)-1; m++) {

        int MaxScore = 0; int Score = 0; int m2 = m; int n = 0;
        int alnLen=0; int TotalnLen=0; int endrow=0; int endcol=0;

        while (m2 < (tSeq.L)-1 && n < (qSeq.L)-1) {

            if (qSeq.getSeqData()[n] == tSeq.getSeqData()[m2]) {
                Score += A[n+1][m2+1];
                alnLen++;
            }

            else {
                if (Score > MaxScore) {
                    MaxScore = Score; Score = 0;
                    TotalnLen = alnLen; alnLen = 0;
                    endrow = n+1; endcol = m2+1;
                }
            }
            m2++; n++;
        }

        MaxScore=Score; TotalnLen=alnLen; endrow=n+1; endcol=m2+1; // necessary to get the last values

        if (MaxScore != 0) {
            maxScorelist.push_back({MaxScore, TotalnLen, endrow, endcol});
        }
    }
    // sort maxScorelist (sort by MaxScore column)
    sort(maxScorelist.rbegin(), maxScorelist.rend());
    // removing duplicates
    maxScorelist.erase(unique(maxScorelist.begin(), maxScorelist.end()), maxScorelist.end());
    // select top 4 ungapped alignments
    vector<vector<int>> top4ungapList(maxScorelist.begin(), maxScorelist.begin() + 4);
    return top4ungapList;
}


struct alignmentInfo completingAlignment (vector<vector<int>> scoreMat, vector<vector<int>> maxScorelist) {

    // finding optimal alignment, starting from the maximum value coord.
    int optScore = 0;
    vector<vector<int>> posVector = {{1, 1},
                                     {1, 0},
                                     {0, 1}};
    vector<int> scorelist{};
    vector<string> cigarlist{"M", "I", "D"};
    string frontcigar;
    string endcigar;

    // going to the upper-left
    unsigned int startrow = (maxScorelist[0][2]) - (maxScorelist[0][1]);
    unsigned int startcol = (maxScorelist[0][3]) - (maxScorelist[0][1]);
    while (startrow > 0 && startcol > 0) {
        scorelist = {scoreMat[startrow-1][startcol-1], scoreMat[startrow-1][startcol],
                     scoreMat[startrow][startcol-1]};
        // finding & marking maximum score position
        int max = 0;
        int maxpos = 0;
        for (size_t a = 0; a < scorelist.size(); a++) {
            if (max < scorelist[a]) {
                max = scorelist[a];
                maxpos = a;
            }
        }
        optScore += max;

        startrow -= posVector[maxpos][0];
        startcol -= posVector[maxpos][1];
        frontcigar += cigarlist[maxpos];
        scorelist.clear();
    }
    reverse(frontcigar.begin(), frontcigar.end());

    // going to the bottom-right
    unsigned int endrow = (maxScorelist[0][2]);
    unsigned int endcol = (maxScorelist[0][3]);
    while (endrow < scoreMat.size() && endcol < scoreMat[0].size()) {
        scorelist = {scoreMat[endrow - 1][endcol - 1], scoreMat[endrow - 1][endcol], scoreMat[endrow][endcol - 1]};
        // finding & marking maximum score position
        int max = 0;
        int maxpos = 0;
        for (size_t a = 0; a < scorelist.size(); a++) {
            if (max < scorelist[a]) {
                max = scorelist[a];
                maxpos = a;
            }
        }
        optScore += max;

        endrow += posVector[maxpos][0];
        endcol += posVector[maxpos][1];
        endcigar += cigarlist[maxpos];
        scorelist.clear();
    }

    //connecting the pieces
    string ungapcigar;
    for (int i=0; i < maxScorelist[0][1]; i ++) {
        ungapcigar += "M";
    }
    string fin_cigar = frontcigar + ungapcigar + endcigar;
    int fin_optScore = optScore + maxScorelist[0][0];
    struct alignmentInfo optAlnResult = {fin_optScore, startrow, endrow, startcol, endcol, fin_cigar};
    return optAlnResult;
}

//function that collects residue index and C-alpha coordinates of the perfectly aligned sequences
void PareunAlign::AlignedResidueIndex(Matcher::result_t  & optAlnResult, int * ires) {
    int qPos = optAlnResult.qStartPos;
    int tPos = optAlnResult.dbStartPos;
    string cigarString = optAlnResult.backtrace;
    for (size_t btPos = 0; btPos < cigarString.size(); btPos++) {
        if (cigarString[btPos] == 'M') {
            ires[qPos] = tPos;
            qPos++;
            tPos++;
        } else if (cigarString[btPos] == 'I') {
            tPos++;
        } else {
            qPos++;
        }
    }
}

string PareunAlign::backtrace2cigar (string backtrace) {
    string cigar; int count=1;
    for (size_t i=1; i < backtrace.size(); i++) {
        if (backtrace.at(i) != backtrace.at(i-1)) {
            cigar += to_string(count) + backtrace.at(i-1);
            count = 1;
        }
        else {
            count += 1;
            if (i == backtrace.size()-1) {
                cigar += to_string(count) + backtrace.at(i-1);
            }
        }
    }
    return cigar;
}

Matcher::result_t PareunAlign::align(Sequence &qSeq, Sequence &tSeq, SubstitutionMatrix *subMat,
                                     EvalueComputation evaluer) {


    // Add substitution matrix (3di) parm
    vector<vector<int>> alignMatrix = initializeScoreMatrix(qSeq, tSeq, * subMat);

    // No need to use 3di here in ungapped align.
    vector<vector<int>> top4ungapFin = ungappedAlign(qSeq, tSeq, alignMatrix);

    // complete the alignment by connecting the ungapped local alignments individually.
    struct alignmentInfo optAlnResult = completingAlignment(alignMatrix, top4ungapFin);

    Matcher::result_t result;
    result.backtrace = optAlnResult.cigar_string;
    result.score = optAlnResult.score;
    result.qStartPos = optAlnResult.query_start;
    result.qEndPos = optAlnResult.query_end;
    result.dbEndPos = optAlnResult.target_end;
    result.dbStartPos = optAlnResult.query_start;
    //result.qCov = SmithWaterman::computeCov(result.qStartPos1, result.qEndPos1, querySeqObj->L);
    result.qcov = SmithWaterman::computeCov(result.qStartPos, result.qEndPos, qSeq.L);
    //result.tCov = SmithWaterman::computeCov(result.dbStartPos1, result.dbEndPos1, targetSeqObj->L);
    result.dbcov = SmithWaterman::computeCov(result.dbStartPos, result.dbEndPos, tSeq.L);
    result.eval = evaluer.computeEvalue(result.score, qSeq.L);

    return result;
}



