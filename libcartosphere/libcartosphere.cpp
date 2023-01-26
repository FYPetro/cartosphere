//
//  libcartosphere.cpp
//  libcartosphere
//
//  Created by Ziqiang Li on 1/25/23.
//

#include <iostream>
#include "libcartosphere.hpp"
#include "libcartospherePriv.hpp"

void libcartosphere::HelloWorld(const char * s)
{
    libcartospherePriv *theObj = new libcartospherePriv;
    theObj->HelloWorldPriv(s);
    delete theObj;
};

void libcartospherePriv::HelloWorldPriv(const char * s) 
{
    std::cout << s << std::endl;
};

