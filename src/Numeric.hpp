//
// Created by will on 9/17/25.
//

#ifndef PDENCLOSE_NUMERIC_H
#define PDENCLOSE_NUMERIC_H
#include <cstdint>
#include <ostream>
#include <vector>

/**
 * Numeric data suitable for approximation schemes.
 */
template<typename T>
concept Numeric = requires(T a, uint32_t nat, double real, std::ostream &out) {
    /*
     * Numeric-numeric operations
     */
    { a + a }       -> std::convertible_to<T>;
    { a - a }       -> std::convertible_to<T>;
    { a / a }       -> std::convertible_to<T>;
    { a * a }       -> std::convertible_to<T>;

    /*
     * Unary operations
     */
    { a.pow(nat) }  -> std::convertible_to<T>;
    { a.abs() }     -> std::convertible_to<T>;
    { a.exp() }     -> std::convertible_to<T>;
    { a.tanh() }     -> std::convertible_to<T>;
    { a.sigmoid() }  -> std::convertible_to<T>;

    /*
     * Comparison operators
     */
    { a == a }      -> std::convertible_to<bool>;
    { a != a }      -> std::convertible_to<bool>;
    { a < a  }      -> std::convertible_to<bool>;
    { a <= a }      -> std::convertible_to<bool>;
    { a > a  }      -> std::convertible_to<bool>;
    { a >= a }      -> std::convertible_to<bool>;

    { a < real  }   -> std::convertible_to<bool>;
    { a <= real }   -> std::convertible_to<bool>;
    { a > real  }   -> std::convertible_to<bool>;
    { a >= real }   -> std::convertible_to<bool>;

    /*
     * Numeric-scalar operations
     */
    { a * real }    -> std::convertible_to<T>;
    { a + real }    -> std::convertible_to<T>;
    { a - real }    -> std::convertible_to<T>;
    { a / real }    -> std::convertible_to<T>;

    /*
     * Abstraction operations
     * Our real-valued wrapper supports aliases for these.
     */
    { a.radius() }  -> std::convertible_to<double>;
    { a.min()    }  -> std::convertible_to<double>;
    { a.max()    }  -> std::convertible_to<double>;

    /*
     * Compositional operations
     * Needed for sound bounds in finite volume methods.
     */
    { a.union_with(a) } -> std::convertible_to<T>;
    { a.split(nat)    } -> std::convertible_to<std::vector<T>>;

    /*
     * Constructors
     */
    { T(real) }         -> std::convertible_to<T>;

    /*
     * Utility operations
     */
    { out << a   }      -> std::convertible_to<std::ostream&>;
};

#endif //PDENCLOSE_NUMERIC_H