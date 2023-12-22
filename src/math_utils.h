/////////////////////////////////////////////////////////////////////////////////////////
#ifndef __MATH_UTILS__
#define __MATH_UTILS__
/////////////////////////////////////////////////////////////////////////////////////////
#include <exception>
#include <string>
/////////////////////////////////////////////////////////////////////////////////////////

class overflow_exception : public std::exception {
private:
    std::string message;
public:
    overflow_exception(const char* message="");
    overflow_exception(const std::string& message);
    virtual const char* what() const throw() override;
};

/////////////////////////////////////////////////////////////////////////////////////////


// static inline int gcd(int a, int b) {
//     assert(a>=0 && b>=0);
//     if (b == 0) {
//         return a;
//     } else {
//         return gcd(b, a % b);
//     }
// }
static inline int gcd(int a, int b) {
    int temp;
    while (b != 0) {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

static inline int lcm(int a, int b) {
    int gcd_val = gcd(a, b);
    int lcm_val = (a * b) / gcd_val;
    return lcm_val;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Exact math operation with overflow check
/////////////////////////////////////////////////////////////////////////////////////////

// throw if addition generates an overflow
inline int add_exact(int a, int b) { 
    long r = long(a) + long(b);
    if (r != int(r))
        throw overflow_exception();
    return r; 
}

// throw if subtraction generates an overflow
inline int subtract_exact(int a, int b) { 
    long r = long(a) - long(b);
    if (r != int(r))
        throw overflow_exception();
    return r; 
}

// throw if multiplication generates an overflow
inline int multiply_exact(int a, int b) { 
    long r = long(a) * long(b);
    if (r != int(r))
        throw overflow_exception();
    return r; 
}

// return -1,0,+1 depending of sign of v
inline int sign3(int v) {
    return (v==0) ? 0 : ((v>0) ? 1 : -1);
}

// return true if a and b are comparable in the less-equal-squared sense
//  - positives and negatives are incomparable
//  - 0 compares with both positives and negatives
//  - 0 compares with 0
inline bool comparable_signs(const int a, const int b) {
    return (sign3(a) * sign3(b)) >= 0;
}

// return true if a <= b in the less-equal-squared sense:
//  - |a| <= |b|   and   a,b signs are comparable 
inline bool less_equal_squared(const int a, const int b) {
    return (abs(a) <= abs(b)) && comparable_signs(a, b);
}

/////////////////////////////////////////////////////////////////////////////////////////
#endif // __MATH_UTILS__