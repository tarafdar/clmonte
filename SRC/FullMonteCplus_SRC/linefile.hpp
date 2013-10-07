/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#ifndef LINEFILE_INCLUDED
#define LINEFILE_INCLUDED
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

class LineFile;

// Provides a crude parser which keeps track of what line and column it's on

class LineFile : protected stringstream {
    int Nl;
    int col;
    int Nerr,Nwarn,Ninfo;
    char commentdelim;
    string fn;
    ifstream is;
    ostream& errstream;

    public:
    using stringstream::ignore;
    enum LFEnumType { LF_EOL=0, LF_EOF=1, LF_SPACES=2, LF_EOLS=3 };
    static const LFEnumType end_line, end_file, skip_space, skip_nls;

    bool eof() const { return is.eof(); }
    bool eol() const { return stringstream::eof(); }

    LineFile(const string& fn_,char commentdelim_='\0',ostream& errstream_=cerr) : Nl(0),col(1),Nerr(0),Nwarn(0),Ninfo(0),commentdelim(commentdelim_),
        fn(fn_),is(fn_.c_str(),ios_base::in),errstream(errstream_)
    {
        if (is.good())
            getNextLine();
        else
            throw string("Failed to open");
    }

    string filePos() const {
        stringstream ss,ss2;
        ss << Nl << '-' << col;
        ss2 << fn << ':' << setw(10) << ss.str() << "  ";
        return ss2.str();
    }

    ostream& error() {
        ++Nerr;
        return errstream << filePos() << "Error   | ";
    }

    ostream& warning() {
        ++Nwarn;
        return errstream << filePos() << "Warning | ";
    }

    ostream& info() {
        ++Ninfo;
        return errstream << filePos() << "Info    | ";
    }

    void processComment();
    void getNextLine();
    void skipWhitespace();

    template<class T>friend LineFile& operator>>(LineFile&,T&);

    friend LineFile& operator>>(LineFile&,LFEnumType);
};

string truncstr(string& s,unsigned maxlen);
LineFile& operator>>(LineFile& lf,LineFile::LFEnumType t);

template<class T>LineFile& operator>>(LineFile& lf,T& x)
{
    stringstream& ss=lf;
    if (ss.eof())
        lf.error() << "Found end of line, expecting more input" << endl;
    else {
        ss >> x;
        if (lf.eof()) { lf.col=ss.tellg(); }
        else if (lf.fail() || (!ss.eof() && !isspace(lf.peek())))
        {
            lf.clear();
            if (lf.tellg() != -1)
                lf.col = lf.tellg();
            string s;
            ss >> s;
            lf.error() << "Invalid input \"" << truncstr(s,20) << "\"" << endl;
        }
        else
            lf.col=ss.tellg();
    }
    return lf;
}
#endif
