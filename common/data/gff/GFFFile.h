/*
 * ============================================================================
 *
 *       Filename:  GFFReader.h
 *
 *    Description:  Define a Naive GFF file for PacBio adapter annotations.
 *
 *        Version:  1.0
 *        Created:  01/09/2014 05:58:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */
#include <string>
#include <iostream>
#include <vector>
using namespace std;

class GFFEntry {
public:
    string name, type, source;
    UInt start, end;
    char strand;
    float score;
    string frame;
    string attributes;
    GFFEntry(string & _name, string & _source, string & _type, 
            UInt & _start, UInt & _end, float & _score, 
            char & _strand, string & _frame, string _attributes) {
        name = _name; 
        source = _source;
        type = _type;
        start = _start;
        end = _end;
        score = _score;
        strand = _strand;
        frame = _frame;
    }
};

class GFFFile {
public:
    vector<GFFEntry> entries;
    void ReadAll(string & gffFileName) {
	    ifstream gffIn;
        CrucialOpen(gffFileName, gffIn, std::ios::in);
        while(gffIn) {
            string line;
            getline(gffIn, line);
            stringstream linestrm(line);
            string name, source, type;
            UInt start, end;
            char strand;
            float score;
            string frame, attributes;
           // A sample record in adapterGffFile:
           // ref000001   .   adapter 10955   10999   0.00    +   .   xxxx
            linestrm >> name >> source >> type 
                     >> start >> end >> score 
                     >> strand >> frame >> attributes;
            entries.push_back(GFFEntry(
                name, source, type, start, end, 
                score, strand, frame, attributes));
        }
        gffIn.close();
    }
};
