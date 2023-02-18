#ifndef math_h
#define math_h


// returns a random number between (including) 0 and 1
double rand01();

// returns a random integer between (including) lower and 
// (excluding) upper bounds.
// Example: If you want to sample the numbers 0, 1, and 2
// randomly, call rand_int(0,3)
int rand_int(int lower_bound, int upper_bound);


#endif
