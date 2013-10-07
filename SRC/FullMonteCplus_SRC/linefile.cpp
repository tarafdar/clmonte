/* This source file is part of FullMonte v0.1
    Copyright (c) Jeffrey Cassidy, 2013

    www.eecg.utoronto.ca/~cassidy/fullmonte

    jeffrey.cassidy@mail.utoronto.ca

    Distributed under the BSD 3-clause license; For details, see LICENSE.txt */

#include "linefile.hpp"

const LineFile::LFEnumType LineFile::skip_nls=LF_EOLS;
const LineFile::LFEnumType LineFile::skip_space = LF_SPACES;
const LineFile::LFEnumType LineFile::end_line = LF_EOL;
const LineFile::LFEnumType LineFile::end_file = LF_EOF;

string truncstr(string& s,unsigned maxlen)
{
    if (s.length() < maxlen)
        return s;
    else
        return s.substr(0,maxlen-3) + "...";
}

void LineFile::skipWhitespace()
{
    while(isspace(peek()) && !eof())
    {
        ++col;
        ignore(1);
    }
}

void LineFile::processComment()
{
    string rest;
    std::getline((stringstream&)*this,rest);
    info() << "Comment: " << rest << endl;
}

void LineFile::getNextLine()
{
    string rest;

    while(!is.eof() && is.peek() == commentdelim)
    {
        ++Nl;
        std::getline(is,rest);
        info() << "Comment: " << rest << endl;
        col=1;
    }
    ++Nl;

    if (!is.eof()){
        is.clear();
        std::getline(is,rest);
        col=1;
        clear();
        str(rest);
    }
}

LineFile& operator>>(LineFile& lf,LineFile::LFEnumType t)
{
    string rest;
    switch(t){
        case LineFile::LF_EOL:
        lf.skipWhitespace();
        if (lf.is.eof()){ lf.error() << "Premature EOF" << endl; }
        else if(lf.eol()){}
        else if(lf.commentdelim && lf.peek() == lf.commentdelim){ lf.processComment(); }
        else
        {
            getline(lf,rest);
            lf.error() << "Expecting newline, found: " << truncstr(rest,20) << endl;
        }
        lf.getNextLine();
        break;

        case LineFile::LF_EOF:
        lf.skipWhitespace();
        if (!lf.is.eof())
            lf.error() << "Expecting end of file" << endl;
        break;

        case LineFile::LF_EOLS:
        lf.skipWhitespace();
        while(lf.is.peek() == '\n')
        {
            ++lf.Nl;
            lf.is.ignore(1);
        }
        lf.getNextLine();
        break;

        case LineFile::LF_SPACES:
        lf.skipWhitespace();
        break;

        default:
        break;
    }
    return lf;
}
