#ifndef DATE_AND_TIME_OF_CURRENT_SYSTEM_STATE_HH
#define DATE_AND_TIME_OF_CURRENT_SYSTEM_STATE_HH

#include <iostream>
#include <sstream> 
#include<fstream>
#include <string>
#include <stdio.h>
#include <ctime>


namespace Adonis{

  //! <a href="http://www.cplusplus.com/reference/ctime/localtime/">Link</a>
  inline std::string show_time(){
    std::time_t rawtime;
    struct tm* timeinfo;
    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    
    std::string str = "Current local time and date: " + std::string(asctime(timeinfo));
    return str;
  }
 

  inline void local_time(){
    std::cout << show_time() << std::endl;
  }

}

#endif
