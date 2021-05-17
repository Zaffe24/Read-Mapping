'''The Burrows-Wheeler transform of a string and string matching'''
'''by Pietro Zafferani'''

import random

'''Creates a random DNA sequence of desired length.'''


def getString(length: int) -> str:
    string = ''.join(random.choice('ACTG') for i in range(length))
    return string


'''Adds the terminal sign to a string.'''


def add_terminal(string: str) -> str:
    return string + '$'


'''It returns a list ordered alphabetically that keeps tracks of the occurrences of each unique element in the string.
    This allows not to store the F column, in fact we use the registry to move along the L column. '''


def occurences(string: str) -> list:
    T = add_terminal(string)
    # create an empty dictionary used to store the appearances of each element
    registry = {}
    for i in list(T):
        if i not in registry.keys():
            registry[i] = 1
        else:
            registry[i] += 1

    # the dictionary is converted in an ordered list
    return sorted([(char, registry[char]) for char in registry], key=lambda j: j[0])


'''Creates the Burrows-Wheeler Transform of the given string.'''


def getBWT(string: str) -> list:
    T = add_terminal(string)
    # create the matrix representing the BWT
    BWT_matrix = []
    for i in range(len(T)):
        # create all the cyclic permutations of the string and append them to the matrix
        # it adds also the offset for each permutation that is needed for the string matching
        BWT_matrix.append([T[i:] + T[:i], i])

    # the rows are ordered alphabetically
    return sorted(BWT_matrix)


'''Prints the BWT properly.'''


def print_BWT(M: list) -> print:
    for row in M:
        print(row)


'''Creates the L column from the BWT and it associates the respective offset.'''


def get_L(M: list) -> list:
    L = []
    for row in M:
        L.append([row[0][-1], row[1]])
    # each element in L is formed by ['character', offset]
    return L


'''It takes as parameter the original string and the occurrences registry, then it adds the indexes to each element in L
in order to simplify the FM Index. '''


def getRankedL(string: str, registry: list) -> list:
    # creates the L column
    L = get_L(getBWT(string))
    # associate an index to each element
    for element in registry:
        index = 0
        for position in L:
            if element[0] == position[0]:
                position.insert(1, index)
                index += 1

    # each rank of L column returned is formed by: ['character', index, offset]
    return L


'''Given a character in L and its index, this function allows to retrieve the successive element in the compressed L
    column. It works by looking in the ordered registry for all the elements that come before the given one and return
    their number, this number correspond to the rank of the next element in L.'''


def get_rank(character: str, index: int, registry: list):
    total = 0
    Range = 0
    # this variable defines whether the elements is present in the registry or not
    present = True
    # iterate over the registry to count all the elements before 'character'
    for element in registry:

        if element[0] < character:
            # the number of previous elements is added to teh count
            total += element[1]
            present = False

        # the given element is found in the register
        elif element[0] == character:
            # this variable tells how many occurrences of the elements there are
            Range = element[1]
            present = True
            break

    if present:
        # the element is found so we return the final rank and its range
        return total + index, Range
    # the element is not present so it's useless to continue the indexing in L
    return 'No matches'


'''It carries out the reverse of the BWT. The only parameter given is the original string.'''


def reverseBWT(string: str) -> str:
    # the reverse string is initiated by default with the terminal sign since the output will be reversed
    reversed = '$'
    occurence = occurences(string)
    L = getRankedL(string, occurence)
    # we start searching the matrix by starting from the first element in L
    element = L[0]
    # the loop is executed until the whole BWT matrix has been searched, namely until the terminal sign is encountered
    while element[0] != '$':
        # add the current element to the reversed string
        reversed += element[0]
        # compute the rank of the next element to be searched
        next_rank = get_rank(element[0], element[1], occurence)
        # update the element with the newly found one
        element = L[next_rank[0]]

    # the string is reversed in order to be like the original one
    return reversed[::-1]


'''Recursive function that aligns the query string to the target starting from a specific rank in L.
    It is based on the concept of FM Index'''


def recursion(start: list, L: list, p: str, registry: list) -> str:
    # base case of recursion when the match is perfect
    if len(p) == 0:
        return 'match'

    # base case of recursion when the match is not perfect
    elif p[-1] != start[0]:
        return 'No match'

    # compute the next element to be checked in L
    rank = get_rank(start[0], start[1], registry)
    # update the parameter with the new rank
    next_el = L[rank[0]]

    # recursive case checks the compatibility of next character in the query by FM Index
    return recursion(next_el, L, p[:-1], registry)


'''Given only the target-string and a query-string, it returns the offsets of all the possible alignments of the 
    query on the target. It is based on the FM Index strategy implemented by a recursive function.'''


def matching_offsets(string: str, p: str):
    registry = occurences(string)
    L = getRankedL(string, registry)
    # we start by defining the range of elements in L that correspond to the last character of the query string
    start = get_rank(p[-1], 0, registry)
    # check whether the last character of the query is present in the registry
    if start != 'No matches':
        # this list will contain the offsets values
        offsets = []
        # we check only the ranks on L that start with the desired character
        for chars in L[start[0]:start[0] + start[1]]:
            # storing the offset value
            offset = chars[2]
            # invoke the recursive function that aligns the query to a sequence in L
            recursiveF = recursion(chars, L, p[:-1], registry)

            # the match is perfect and thus it can be saved
            if recursiveF == 'match':
                offsets.append(offset)

        return sorted(offsets, reverse=False)

    # eventuality that the query is not present in the target
    return None


'''Prints all the possible alignments, if there are any, of the query on the target string.'''


def align_matches(string: str, p: str, offsets) -> print:
    # check whether any offsets are presents
    if offsets != None:

        # print each alignment and its relative range on the target string
        for n in offsets:
            space = (n - len(p) + 1)
            print(str('range-> ' + '[' + str(space) + ':' + str(space + len(p)) + ']'))
            print(' ' * space + p)
            print(string + '\n')
    else:
        # no offsets found means that the query is not present in the target
        print('read/s not mappable')


'''Function needed for testing the code as required by the assignment.'''


def insertRandom(string: str, p: str, n: int) -> str:
    lp = len(p)
    ls = len(string)
    # depends on how many times we want to perform the insertion
    for i in range(n):
        # create a random offset
        Rindex = random.randrange(0, ls - 1, lp)
        # insert the query string at the offset obtained
        string = string[:Rindex] + p + string[Rindex:]

    # the modified string is returned
    return string


'''First testing function: takes one query string to insert in the target and then returns the mapping.
    Mapping is always possible'''


def Test1():
    lenS = input('Please choose the genome\'s length: ')
    String = getString(int(lenS))
    print(String + '\n')
    P = str(input('Please type in the query sequence to insert in the genome: ')).upper()
    newS = insertRandom(String, P, 1)
    print(newS + '\n')
    Matchings = matching_offsets(newS, P)
    align_matches(newS, P, Matchings)


'''Second testing function: takes one query string to insert 2 times in the target and then returns the mapping.
    Mapping is always possible'''


def Test2():
    lenS = input('Please choose the genome\'s length: ')
    String = getString(int(lenS))
    print(String + '\n')
    P = str(input('Please type in the query sequence that will be inserted 2 times in the genome: ')).upper()
    newS = insertRandom(String, P, 2)
    print(newS + '\n')
    Matchings = matching_offsets(newS, P)
    align_matches(newS, P, Matchings)


'''Third testing function: creates a random target string, the user is asked to map a query string that must not be
    contained in the target. Obviously, the mapping's outcome should be negative in this case.'''


def Test3():
    lenS = input('Please choose the genome\'s length: ')
    String = getString(int(lenS))
    print(String + '\n')
    P = str(input('Please type in a query sequence that is NOT present in the genome: ')).upper()
    Matchings = matching_offsets(String, P)
    print()
    align_matches(String, P, Matchings)


'''It is the function that interacts with the user and allows to choose what test to perform.'''


def MAIN() -> print:

    I = eval(input('Choose what type of test would you like to perform: [1,2,3]? '))

    if I == 1:
        return Test1()  # map one read to genome
    elif I == 2:
        return Test2()  # two identical reads to genome
    else:
        return Test3()  # the mapping is not possible


if __name__ == '__main__':
    # run the file to start the testing
    MAIN()



    # print_BWT(getBWT(string)) -> prints the BWT matrix

    # reverseBWT(string) -> returns the original string from the BWT reversal

    # matching_offfsets(string, p) -> returns a list containing the positions in which the reads align to the target

    # align_matches(string, p, offsets) -> prints the correct mapping of the query on the target together with the
    # respective numerical range
