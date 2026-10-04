#! /usr/bin/env bash
# Brad Bell's reduced test, coin-or/CppAD#259 comment 5960186551
# (attachment swap.sh). Body below is unchanged.
# https://github.com/coin-or/CppAD/issues/259#issuecomment-5960186551
set -e -u
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# -----------------------------------------------------------------------------
#
# temp.cpp
cat << EOF > temp.cpp
#include <utility>
#include <iostream>
//
// myclass
class myclass {
public:
    int m_value;
    myclass(void)
    { }
    myclass(int value) : m_value(value)
    { }
    void swap(myclass&& other)
    {   std::swap(m_value, other.m_value); }
};
//
// add
myclass add(const myclass& left, const myclass right) {
    myclass result;
    result.swap( myclass( left.m_value + right.m_value ) );
    return result;
}
//
// main
int main(void) {
    //
    // left, right
    myclass left(2), right(3);
    //
    // result
    myclass result = add(left, right);
    //
    std::cout 
        << left.m_value << "+" 
        << right.m_value << " = " 
        << result.m_value << "\n";
    //
    return 0;
}
EOF
#
# temp
echo_eval g++ temp.cpp -o temp \
    -O2 \
    -Wall \
    -Wextra \
    -Wpedantic \
    -Wshadow \
    -Wconversion \
    -Wlogical-op \
    -Wduplicated-cond \
    -Wduplicated-branches \
    -Wunused \
    -Wold-style-cast \
    -Woverloaded-virtual \
    -Wnull-dereference \
    -Wformat=2 -Werror
#
echo_eval ./temp

