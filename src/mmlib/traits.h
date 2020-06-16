#ifndef __TRAITS_H
#define __TRAITS_H

/// trait used in nonbonded_periodic.h to select which block 
/// of inline code to compile in
template< typename T_Type1, typename T_Type2 >
struct is_same_type{ static const bool value = false; };
/// This template specialisation essentially acts as a compile-time
/// if statment: if you access is_same_type<typea,typeb>::value it will
/// be true if typea and typeb are the same Type and false otherwise.
/// used in an if statement the compiler can compile away the if
/// because the result is already known at compile time (as if you
/// said if(true){ ... }
template< typename T >
struct is_same_type<T,T>{
	static const bool value = true;
}; 

#endif

