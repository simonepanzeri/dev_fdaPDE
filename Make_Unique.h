//
// Created by simonepanzeri on 26/01/2022.
//

#ifndef SPACE_TIME_CPP_MAKE_UNIQUE_H
#define SPACE_TIME_CPP_MAKE_UNIQUE_H

#include <memory>

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique_time(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

#endif //SPACE_TIME_CPP_MAKE_UNIQUE_H
