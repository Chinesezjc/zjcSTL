#include "bignumber.hpp"
int main()
{
	zjcSTL::bigint<10> a, b;
	std::cin >> a >> b;
	std::cout << a + b << std::endl;
	return 0;
}