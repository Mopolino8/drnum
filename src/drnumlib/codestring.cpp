#include "codestring.h"


CodeString::CodeString()
{
  string::operator=("");
}


CodeString::CodeString(const string& inputstring)
{
  istringstream iss_input(inputstring);
  buildFrom(iss_input);
}


CodeString::CodeString(istringstream& iss_input)
{
  buildFrom(iss_input);
}


void CodeString::buildFrom(istringstream& iss_input)
{
  // Read input istringstream until its end and copy to string in "this"
  //  * treat '/n' ',' ';' as delimiter ' '
  //  * allow only one ' ' as delimiter
  //  * eliminate leading and trailing blancs
  //  * eliminate leading zeroes in code string and let zero code be 1-digit
  // Example:
  // input "  aa, 034  00\n 70;5 " converts to "aa 34 0 70 5"

  string::operator=(""); // empty
  int blanccount = 0;
  while (iss_input.good())  // loop while extraction from file is possible
  {
    char c = iss_input.get();
    if(c == '\n' || c == ',' || c == ';') { // treat as delimiter ' '
      c = ' ';
    }
    if(c == ' ') { // allow only one ' ' as delimiter
      blanccount++;
    } else {
      blanccount = 0;
    }
    if(blanccount < 2) {
      if (iss_input.good()) { // avoid end marker
        push_back(c);
      }
    }
  }
  // eliminate leading and trailing blancs, if applicable
  if(*(begin()) == ' ') {
    erase(begin());
  }
  size_t len = length();
  if(len > 0) {
    char lastchar = operator[](len-1);
    if(lastchar == ' ') {
      erase(len-1);
    }
  }
  // eliminate leading zeroes in code words
  bool startacode = true;
  size_t is = 0;
  while(is < length()-1) {
    char c0 = operator[](is);
    char c1 = operator[](is+1);
    if(c0 == '0' && c1 != ' ' && startacode) {
      erase(is,1);
    }
    else if(c0 == ' ') {
      startacode = true; // a new code starts in next iteration
      is++;
    }
    else {
      startacode = false; // inside character sequence of a code
      is++;
    }
  }
}


CodeString CodeString::compareCodeString(const CodeString& otherstring,
                                         bool& allsame,
                                         bool& overruled,
                                         bool& underruled)
{
  allsame = false;
  underruled = false;
  overruled = false;
  CodeString concatcodestring;

  // check, if strings these are fully equal
  if(compare(otherstring) == 0) {
    allsame = true;
    overruled = false;
    underruled = false;
    concatcodestring = otherstring;
    return concatcodestring;
  }

  // check overruled / underruled
  bool totaldiff = false;
  //.. transfer to istringstreams to facilitate
  istringstream iss_this(*this);
  istringstream iss_other(otherstring);
  //.. check word by word
  while(iss_this.good() && iss_other.good()) {
    //.... get next code words from both streams
    string word_this, word_other, word_concat;
    getline (iss_this, word_this, ' ');
    getline (iss_other, word_other, ' ');
    bool thisgood = iss_this.good();
    bool othergood = iss_other.good();
    //.... check valid inputs
    if((thisgood && !othergood) || (!thisgood && othergood)) {
      totaldiff = true; // due to different word counts
      break;
    }
    //.... check, if corresponding code words match
    if(word_this.compare(word_other) != 0) {
      if(word_this.compare("0") == 0) {   // check overruled
        overruled = true;
        word_concat = word_other;
      }
      else if (word_other.compare("0") == 0) { // check underruled
        underruled = true;
        word_concat = word_this;
      }
      else {                        // totally different, skip
        totaldiff = true; // due to different code words
        break;
      }
    }
    else {
      word_concat = word_this;
    }
    concatcodestring += ' ';
    concatcodestring += word_concat;
  }
  if(*(concatcodestring.begin()) == ' ') {
    concatcodestring.erase(concatcodestring.begin());
  }
  // handle totally non matching code strings
  if(totaldiff) {
    allsame = false;
    underruled = false;
    overruled = false;
    concatcodestring.string::operator=("");  // invalid anyway
  }
  return concatcodestring;
}


size_t CodeString::countNumCodeWords()
{
  size_t num_words = 0;
  istringstream iss_this(*this);
  // Read word by word and count
  while(iss_this.good()) {
    // read next code word
    string word;
    getline (iss_this, word, ' ');
    num_words++;
  }
  return num_words;
}


string CodeString::accessCodeWord(const size_t& i_word,
                                  bool& success)
{
  success = false;
  string word;
  size_t count_words = 0;
  istringstream iss_this(*this);
  // Read word by word and count
  while(iss_this.good() && count_words <= i_word) {
    // read next code word
    getline (iss_this, word, ' ');
    if(count_words == i_word) {
      success = true;
      break;
    }
    count_words++;
  }
  if(!success) {
    word = "";
  }
  return word;
}


CodeString::~CodeString()
{
  // no own data
}
