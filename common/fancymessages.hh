#ifndef FANCY_MESSAGES_THAT_CAN_BE_LAUNCHED_TO_SCREEN_HH
#define FANCY_MESSAGES_THAT_CAN_BE_LAUNCHED_TO_SCREEN_HH

#include <iostream>
#include <stdio.h>   //for printf
#include <string>
#include <fstream>

#include "date.hh"

namespace Adonis{
  
  // typedef enum{
    
  // }ForegroundBackgroundColors;

  /**
   *  \brief Some cool screen outputs ;-)
   *
   *  \code
      Adonis::FancyMessages eowin; 
      eowin.nice_error();  //use default arguments
      //changes foreground color
      eowin.nice_error("Your error description stands here", eowin.green_FG()); 
   *  \endcode
   */
  class FancyMessages{
  private:
    //!The colors - colors beginning with an upper case letter are FG colors
    enum colors_{orig = 0, Black = 30, Red = 31, Green = 32, Yellow = 33, Blue = 34, Purple = 35, Cyan = 36, White = 37,  
		 black = 40, red = 41, green = 42, yellow = 43, blue = 44, purple = 45, cyan = 46, white = 47};

  public:
    FancyMessages(){}

    /**
     *  \brief Here you can define some colors -- just in case if you, e.g., want to alter the foreground colour but don't know the encodings for them.
     * I just wrote functions for the coolest colors
     */
    inline const int red_FG() { return Red;}
    inline const int green_FG() {return colors_(32);}
    inline const int blue_FG() {return Blue;}
    inline const int purple_FG() {return Purple;} 

    //!the same for the background colours
    inline const int black_BG(){return black;}
    inline const int green_BG(){return green;}
    inline const int white_BG(){return white;}


    /**
     * \brief specify FG/BG colours and message string
     * \param fg Foreground colour as stated in colors_ {default: Red}
     * \param bg Backgroun colour as stated in colors_ {default: black}
     * \param s Message {default: std::string()}
     */

    inline void nice_error(const std::string& s = std::string(),const int fg = Red, const int bg = black ) const{
      
      //!E.g.: 'printf("%c[%d;%dmFuck", 0x1B, 1,   31)' prints red, bright, bold and transparent background
      printf("%c[%d;%d;%dmFuck", 0x1B, 1,  fg,  bg);
    //                         |     |
    //         prepare all colors  text appearance(1 = brighter than 0 and bold
    //                                              4 = dim and underlined)    
    
      std::cout<<"\n"<<std::endl;
      std::cout<<"********  ERROR **** ERROR **** ERROR **** ERROR **** ERROR ******* :-(  "<<std::endl;
	 
    std::cout<<"           ________ "<<std::endl;
    std::cout<<"          /\\  _____\\  "<<std::endl;	
    std::cout<<"          \\ \\ \\____/ "<<std::endl;  
    std::cout<<"           \\ \\ \\______   ____  ____  _____  ____"<<std::endl;
    std::cout<<"            \\ \\  _____\\ /  __`/  __`/  __ `/  __`."<<std::endl;
    std::cout<<"             \\ \\ \\____/ /\\ \\__/\\ \\__/' \\ \\ \\\\ \\__/"<<std::endl;
    std::cout<<"              \\ \\ \\_____\\ \\ \\ \\ \\ \\ \\ \\ \\_\\ \\\\ \\  "<<std::endl;
    std::cout<<"               \\ \\_______\\ \\_\\ \\ \\_\\ \\ \\____/ \\_\\ "<<std::endl;
    std::cout<<"                \\/_______/\\/_/  \\/_/  \\/____/\\/_/ "<<std::endl;

    std::cout<<std::endl;

   

    //print to screen 
  std::cout<< s <<std::endl<<std::endl<<"I STOP ALL COMPUTATIONS, MARC! \n \n";
  
  
  std::cout<<"********  ERROR **** ERROR **** ERROR **** ERROR **** ERROR ******* :-( "<<std::endl;
  //std::cout<<"\n"<<std::endl; 


  
  printf("%c[%dm", 0x1B, orig); //RESET to original console appearance 
    
  
    }
  
    /**
     * \brief Same error function except that the letter 'r' is more prominent!
     * \param fg Foreground colour as stated in colors_ {default: Red}
     * \param bg Backgroun colour as stated in colors_ {default: black}
     * \param s Message {default: std::string()}
     */
    inline void nice_eRRoR(const std::string& s = std::string(),const int fg = Red, const int bg = black) const{
      
      printf("%c[%d;%d;%dmFuck", 0x1B, 1,  fg,  bg);
    //                         |     |
    //         prepare all colors  text appearance(1 = brighter than 0 and bold
    //                                              4 = dim and underlined)    
    
      std::cout<<"\n"<<std::endl;
      std::cout<<"********  ERROR **** ERROR **** ERROR **** ERROR **** ERROR ******* :-(  ";
	 
    //cool output ;-)
    //NOTE: a backslash ('\') is given via '\\' since the first one is a reserved character!
    std::cout<< std::endl;
    std::cout<<"              ________" << std::endl;
    std::cout<<"             /\\  _____\\ " <<std::endl;		
    std::cout<<"             \\ \\ \\____/ " <<std::endl;  
    std::cout<<"              \\ \\ \\______   _____  _____  _____  _____ "<<std::endl;
    std::cout<<"               \\ \\  _____\\ /  __ \\/  __ \\/  __ `/  __ \\ "<<std::endl;
    std::cout<<"                \\ \\ \\____/ /\\ \\/\\_\\\\ \\/\\_\\' \\ \\ \\\\ \\/\\_\\ "<<std::endl;
    std::cout<<"                 \\ \\ \\_____\\ \\ \\/_/ \\ \\/_/ \\ \\_\\ \\\\ \\/_/ "<<std::endl;
    std::cout<<"                  \\ \\_______\\ \\_\\  \\ \\_\\  \\ \\____/ \\_\\ "<<std::endl;
    std::cout<<"                   \\/_______/\\/_/   \\/_/   \\/____/\\/_/ "<<std::endl;
    
    std::cout<<std::endl;  
   

    //print to screen 
  std::cout<< s <<std::endl<<std::endl<<"I STOP ALL COMPUTATIONS, MARC! \n \n";
  
  
  std::cout<<"********  ERROR **** ERROR **** ERROR **** ERROR **** ERROR ******* :-( "<<std::endl;
  //std::cout<<"\n"<<std::endl; 


  
  printf("%c[%dm", 0x1B, orig); //RESET to original console appearance 
    
  
    }



    /**
     * \brief Outputs an information embedded in a nice setting, which can be used, e.g. when some threshold is overshot, some iteration counter is exceeded, et cetera. A foretaste is given below
     
           _______           _____
          /\__  __\         /  ___`.
          \/_/\ \_/  __     /\ \___/
             \ \ \  /\ \____\ \ \__  ____ 
              \ \ \ \ \  __ `\ \  _\· __ `.
               \ \ \ \ \ \_/\ \ \ \/\ \ \ \
                \_\ \_\ \ \\ \ \ \ \ \ \_\ \  
               /\ _____\ \_\\ \_\ \_\ \ ___/  
               \/______/\/_/ \/_/\/_/\/____/

     
     */
    inline void nice_info(const std::string& s = std::string(),const int fg = Blue, const int bg = black) const{
      
      printf("%c[%d;%d;%dmAttention", 0x1B, 1,  fg,  bg);
    //                         |     |
    //         prepare all colors  text appearance(1 = brighter than 0 and bold
    //                                              4 = dim and underlined)    
    
      std::cout<<"\n"<<std::endl;
      std::cout<<"********  INFO **** INFO **** INFO **** INFO **** INFO **** INFO ******* :-/ ";
	 
    //cool output ;-)
    //NOTE: a backslash ('\') is given via '\\' since the first one is a reserved character!
    std::cout<< std::endl;
    std::cout<<"              _______           _____ "<<std::endl;
    std::cout<<"             /\\__  __\\         /  ___`. "<<std::endl;
    std::cout<<"             \\/_/\\ \\_/  __     /\\ \\___/ "<<std::endl;
    std::cout<<"                \\ \\ \\  /\\ \\____\\ \\ \\__  ____ "<<std::endl;
    std::cout<<"                 \\ \\ \\ \\ \\  __ `\\ \\  _\\· __ `. "<<std::endl;
    std::cout<<"                  \\ \\ \\ \\ \\ \\_/\\ \\ \\ \\/\\ \\ \\ \\ "<<std::endl;
    std::cout<<"                   \\_\\ \\_\\ \\ \\\\ \\ \\ \\ \\ \\ \\_\\ \\ "<<std::endl;  
    std::cout<<"                  /\\ _____\\ \\_\\\\ \\_\\ \\_\\ \\ ___/ "<<std::endl;  
    std::cout<<"                  \\/______/\\/_/ \\/_/\\/_/\\/____/ "<<std::endl;

    
    std::cout<<std::endl;  
   

    //print to screen 
  std::cout<< s <<std::endl<<std::endl;
  
  
  std::cout<<"********  INFO **** INFO **** INFO **** INFO **** INFO **** INFO ******* :-/ "<<std::endl;
  //std::cout<<"\n"<<std::endl; 


  
  printf("%c[%dm", 0x1B, orig); //RESET to original console appearance 
    
  
    } 
    



    /**
     * \brief Outputs a warning embedded in a nice setting, which can be used, e.g. when a numerical calculation succeeded but yields no meaningful values, etc.
     
   __       __ 
  /\ \     /\ \     
  \ \ \    \ \ \                   __      __  __ 
   \ \ \    \ \ \    ____.   _____/\ \____/\_\/\ \____   _____
    \ \ \    \_\ \  / ___ \ /   __`.\  __ `\_/\ \  __ `\/ ___ \
     \ \ \  / . \ \/\ \  \ \/\ \__/\ \ \_/\ \\`\ \ \_/\ \ \  \ \   
      \ \ \/ / \ \ \ \ \__\ \ \ \   \ \ \\ \ \\ \ \ \\ \ \ \__\ \
       \ \_ / \ \__/\ \____. \ \_\   \ \_\\ \_\\_\ \_\\ \_\ '___ \
        \/_/   \/_/  \/___ /_/\/_/    \/_/ \/_//_/\/_/ \/_/\/___\ \
                                                            /\_____\
                                                            \/_____/  

     
     */
    inline void nice_warning(const std::string& s = std::string(),const int fg = Yellow, const int bg = black) const{
      
      printf("%c[%d;%d;%dmAttention", 0x1B, 1,  fg,  bg);
    //                         |     |
    //         prepare all colors  text appearance(1 = brighter than 0 and bold
    //                                              4 = dim and underlined)    
    
      std::cout<<"\n"<<std::endl;
      std::cout<<"********  WARNING **** WARNING **** WARNING **** WARNING **** WARNING ****** :-| ";
	 
    //cool output ;-)
    //NOTE: a backslash ('\') is given via '\\' since the first one is a reserved character!
      std::cout << "    __       __ " << std::endl; 
      std::cout << "    /\\ \\     /\\ \\" << std::endl;     
      std::cout << "    \\ \\ \\    \\ \\ \\                   __      __  __ "<<std::endl; 
      std::cout << "     \\ \\ \\    \\ \\ \\    ____.   _____/\\ \\____/\\_\\/\\ \\____   _____ " << std::endl;
      std::cout <<"      \\ \\ \\    \\_\\ \\  / ___ \\ /   __`.\\  __ `\\_/\\ \\  __ `\\/ ___ \\ " << std::endl; 
      std::cout << "       \\ \\ \\  / . \\ \\/\\ \\  \\ \\/\\ \\__/\\ \\ \\_/\\ \\\\`\\ \\ \\_/\\ \\ \\  \\ \\" << std::endl;   
      std::cout << "        \\ \\ \\/ / \\ \\ \\ \\ \\__\\ \\ \\ \\   \\ \\ \\\\ \\ \\\\ \\ \\ \\\\ \\ \\ \\__\\ \\ " << std::endl;
      std::cout <<"         \\ \\_ / \\ \\__/\\ \\____. \\ \\_\\   \\ \\_\\\\ \\_\\\\_\\ \\_\\\\ \\_\\ '___ \\ " << std::endl;
      std::cout << "          \\/_/   \\/_/  \\/___ /_/\\/_/    \\/_/ \\/_//_/\\/_/ \\/_/\\/___\\ \\ " << std::endl;
      std::cout <<"                                                              /\\_____\\ " << std::endl;
      std::cout <<"                                                              \\/_____/ " << std::endl;  

     

     
     std::cout<<std::endl;  

    //print to screen 
  std::cout<< s <<std::endl<<std::endl;
  
  
  std::cout<<"********  WARNING **** WARNING **** WARNING **** WARNING **** WARNING ****** :-| ";
  //std::cout<<"\n"<<std::endl; 


  
  printf("%c[%dm", 0x1B, orig); //RESET to original console appearance 
    
  
    } 
    

    
    /**
     * \brief Outputs a warning embedded in a nice setting, which can be used, e.g. when a numerical calculation succeeded but yields no meaningful values, etc.
     
     ________                     __    
    / ______ \                   /\ \
   /\ \   \ \_\                  \_\ \__ __          __
   \ \ \   \/_/     ____. __   __/\__  _\\_\   ____ /\ \____ 
    \ \ \      __  / ___ \\ \ /\ \/_/\ \//_/ /· __ `. \  __ `\
     \ \ \    /\ \/\ \_/\ \\ \\ \ \ \ \ \/\`\/\ \ \ \\ \ \_/\ \
      \ \ \___\_\ \ \ \\_\ \\ \\_\ \ \ \ \ \ \ \ \_\ \\ \ \\ \ \
       \ \________/\ \____. \\_____/  \ \_\ \_\ \ ___/ \ \_\\ \_\
        \/_______/  \/___ /_//____/    \/_/\/_/\/____/  \/_/ \/_/

  
     
     */
    inline void nice_caution(const std::string& s = std::string(),const int fg = Yellow, const int bg = black) const{
      
      printf("%c[%d;%d;%dmAttention", 0x1B, 1,  fg,  bg);
    //                         |     |
    //         prepare all colors  text appearance(1 = brighter than 0 and bold
    //                                              4 = dim and underlined)    
    
      std::cout<<"\n"<<std::endl;
      std::cout<<"********  CAUTION **** CAUTION **** CAUTION **** CAUTION **** CAUTION ****** :-| ";
	 
    //cool output ;-)
    //NOTE: a backslash ('\') is given via '\\' since the first one is a reserved character!
     

      std::cout<<"     ________                     __ "<<std::endl;   
      std::cout << "     / ______ \\                   /\\ \\"<<std::endl;
      std::cout << "    /\\ \\   \\ \\_\\                  \\_\\ \\__ __          __"<<std::endl;
      std::cout << "    \\ \\ \\   \\/_/     ____. __   __/\\__  _\\\\_\\   ____ /\\ \\____ "<<std::endl;
      std::cout<<"     \\ \\ \\      __  / ___ \\\\ \\ /\\ \\/_/\\ \\//_/ /· __ `. \\  __ `\\"<<std::endl;
      std::cout << "      \\ \\ \\    /\\ \\/\\ \\_/\\ \\\\ \\\\ \\ \\ \\ \\ \\/\\`\\/\\ \\ \\ \\\\ \\ \\_/\\ \\"<<std::endl;
      std::cout << "       \\ \\ \\___\\_\\ \\ \\ \\\\_\\ \\\\ \\\\_\\ \\ \\ \\ \\ \\ \\ \\ \\_\\ \\\\ \\ \\\\ \\ \\"<<std::endl;
      std::cout << "        \\ \\________/\\ \\____. \\\\_____/  \\ \\_\\ \\_\\ \\ ___/ \\ \\_\\\\ \\_\\"<<std::endl;
      std::cout << "         \\/_______/  \\/___ /_//____/    \\/_/\\/_/\\/____/  \\/_/ \\/_/"<<std::endl;

     
     std::cout<<std::endl;  

    //print to screen 
  std::cout<< s <<std::endl<<std::endl;
  
  
  std::cout<<"********  CAUTION **** CAUTION **** CAUTION **** CAUTION **** CAUTION ****** :-| ";
  //std::cout<<"\n"<<std::endl; 


  
  printf("%c[%dm", 0x1B, orig); //RESET to original console appearance 
    
  
    } 
    

    /**
     * \brief Just generate colorful output for a std::string
     */
    void nice_output(const std::string& str, const int fg = Green, const int bg = black) const{
      printf("%c[%d;%d;%dm", 0x1B, 1,  fg,  bg);
      std::cout << str << std::endl;
      printf("%c[%dm", 0x1B, orig);  //reset afterwards
    }

    
    //! usual output -- convenient member
    void display(const std::string& str){
      std::cout << str << std::endl;
    }


    void display_in_bold(const std::string& str){
      //!Bash: echo -e "\033[1mtest\033[0m again"
      printf("%c[1m",0x1B);
      std::cout << str << std::endl;
      printf("%c[0m",0x1B);   //reset
    }

    void error_log(const std::string s){
      std::ofstream of("ERROR.log");
      of << "***** automatically generated error file *****"<<std::endl<<std::endl;
      of << s;
      of << std::endl<< std::endl<< show_time() << std::endl;
      of << "**********************************************"<<std::endl;
      of.close();
    }
  };


} //end Adonis

#endif
