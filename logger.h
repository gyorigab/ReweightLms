#ifndef LOGGER_H
#define LOGGER_H

#include <ios>
#include <time.h>
#include <iomanip>

using namespace std;

class Logger
{
public:

    Logger(){}

    enum TYPE{ERROR,WARNING,INFO,DEBUG1,DEBUG2,DEBUG3};

    static void comment(TYPE t,const string& message)
    {
        switch(t)
        {
        case ERROR:                  print(" ERROR    ",message)  ; break;
        case WARNING: if( verb > 0 ) print(" WARNING  ",message)  ; break;
        case INFO:    if( verb > 1 ) print(" INFO     ",message)  ; break;
        case DEBUG1:  if( verb > 2 ) print(" DEBUG1   ",message)  ; break;
        case DEBUG2:  if( verb > 3 ) print(" DEBUG2   ",message)  ; break;
        case DEBUG3:  if( verb > 4 ) print(" DEBUG3   ",message)  ; break;
        default: print(" ",message)  ; break;
        }
    }

    static void print(const string& type, const string& message)
    {
        static const char mon[][4] = {
          "Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
        };

        time_t now = time(0);
        tm localtm = *localtime(&now);

        out << "[LOG] ";

        out << mon[localtm.tm_mon]    << " "
            << localtm.tm_mday        << " "
            << 1900 + localtm.tm_year << " "   ;

        out << setfill('0')
            << setw(2) << localtm.tm_hour << ":"
            << setw(2) << localtm.tm_min  << ":"
            << setw(2) << localtm.tm_sec  << setfill(' ') ;

        out << type << message << "\n";
    }

    void set_verbosity(int fverb){verb = fverb;}

private:

    static int verb;
    static ostream& out;
};

#endif // LOGGER_H
