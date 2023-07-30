/*
	PROTECTED_VECTOR.HPP
	--------------------
*/
#pragma once

#include <stdint.h>

#include <atomic>
#include <vector>


/*
	CLASS PROTECTED_VECTOR
	----------------------
*/
template <typename T>
class protected_vector : public std::vector<T>
	{
	private:
		std::atomic<uint32_t> lock;

	public:
		/*
			PROTECTED_VECTOR::PROTECTED_VECTOR()
			------------------------------------
		*/
		protected_vector() :
			lock(0)
			{
			/* Nothing */
			}

		/*
			PROTECTED_VECTOR::PUSH_BACK()
			-----------------------------
		*/
		void push_back(const T &val)
			{
			uint32_t expected = 0;
//			uint32_t collision = 0;

			while (!lock.compare_exchange_strong(expected, 1))
				{
				expected = 0;
//				collision = 1;
				}

			std::vector<T>::push_back(val);

			lock = 0;
//			if (collision)
//				puts("C");
			}
	};
