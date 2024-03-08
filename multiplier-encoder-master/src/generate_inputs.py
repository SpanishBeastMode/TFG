from itertools import product

def generate_a_prime_number(size):
    lower, upper = 2**(size-1), 2**size-1
    for num in range(upper, lower, -1):
        if num > 1:
            for i in range(2, num):
                if (num % i) == 0:
                    break
            else:
                yield num

# def take_a_few_prime_numbers(size, counter=10):
#     i = 0
#     for num in generate_a_prime_number(size):
#         if i >= counter: break

#         yield num
#         i += 1

def generate_inputs(length, width, counter=1):
    A = generate_a_prime_number(length)
    B = generate_a_prime_number(width)
    P = []
    for a, b in product(A, B):
        out = a*b
        print(f'{a}*{b}={out}')
        P.append(out)
        if len(set(P)) >= counter:
            return set(P)
    return set(P)


if __name__ == '__main__':
    length, width = 3, 4
    inputs = generate_inputs(length, width, counter=1)
    print(inputs)