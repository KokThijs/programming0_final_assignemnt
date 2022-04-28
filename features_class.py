'''
Name: features_class
This module contains my Features class to use in my main module (final_assessment.py).

Author: Thijs Kok
Student number: 348224
Class: Data Science for Life Science

This module contains the class that contains each individual feature and does some work on it.
'''

import re
from genbankparser_class import GenBankParser



class Features(GenBankParser):
    """
    This class gets a list containing the isolated feature lines.

    - function: - feature_to_sequence_builder > this function selects which actions are required
                  to obtain the right sequence string.

                - complement_string > this function returns the reversecomplement of the required
                  sequence

                - join_string > this function seeks what parts are to be joined and does so.

                - return_sequence > the most basic operation. Looks which to sequence indices are
                  required for the sequence associated with the feature and returns a string cont
                  -aining the sequence.

                - complement_join > return the reverse complement of a joined sequence.

                - order_amino_acids > returns an alphabetically ordered string of amino acids

                - sequence_source_isolator > isolated the / 'gene=<insert gene> and / source= ...
                  part of the feature
    """
    mode_select = False # False = seperated, True = uppercased
    object_count = 0 # counts the amount of feature objects created. each init increases count by 1
    features = '' #
    feature_name = ''   # store the name of the feature in each instance
    sequence_source = ''   # dna/rna/protein -sequence
    output_string = '1' # string to be written to the output file


    def __init__(self, feature_string):
        self.features = feature_string.split('/')
        Features.object_count += 1
        self.feature_name = self.name_isolator()
        self.sequence_source = self.sequence_source_isolator()
        self.output_sequence = self.feature_to_sequence_builder()
        self.output_string = self.finished_output()
        self.mode_select = Features.mode_select

    def feature_to_sequence_builder(self):
        """
        Splits the sequence indices and return the corresponding parts of the sequence.

        Input ->    mode is the mode in which to return each string. 1 = seperated, 2 = uppercase
                    standard mode is seperated (i.e. no mode is given)

        Output ->   sequence string in the way it should be expected depending on the mode
        """
        split_first_line = self.features[0].split()
        # joins the seperate entries of the feature indices together to work on them, starts at
        # the second element of the list
        sequence_indices = ''.join(split_first_line[1:])

        # seeks a join statement
        if sequence_indices.startswith('join'):
            return self.string_length_modifier(self.join_string(sequence_indices))

        # seeks if a string contains a join complement AND a join statement
        elif sequence_indices.startswith('complement') and 'join' in sequence_indices:
            # removes the 'complement' and the brackets from the string and saves it
            sequence_indices = sequence_indices.replace('complement(', '').replace(')', '')
            # runs the complement join function which returns the complement of the sequence
            temp_string = self.complement_join(self.join_string(sequence_indices))
            # returns the string in a maximum of 60 characters per line
            return self.string_length_modifier(temp_string)

        # seeks a complement statement
        elif sequence_indices.startswith('complement'):
            # feed the complement_string function with the feat_indices string
            return self.string_length_modifier(self.complement_string(sequence_indices))

        # seeks a string starting with digits
        elif re.search(r'^\d', sequence_indices):
            # returns the output of the return sequence in maxumum of 60 characters per line
            return self.string_length_modifier(self.return_sequence(sequence_indices))

        # checks for an order statement in the protein file
        elif sequence_indices.startswith('order'):
            # remove 'order(' and ')'
            sequence_indices = sequence_indices.replace('order(', '').replace(')', '')
            # joins the indices together and assigns them to temp_string
            temp_string = self.join_string(sequence_indices)
            # orders the string returned by the join function
            temp_string = self.order_amino_acids(temp_string)
            # returns the string with a maximum length of 60 characters
            return self.string_length_modifier(temp_string)

        else:
            print(sequence_indices)
            print('this is not suppposed to be here :(')


    def complement_string(self, index_pair_string):
        '''
        This function takes an index pair string or a sequence string and return the complement
        or reverse complement of the DNA sequence related to it.
        Checks for any > or < markup and handles them accordingly.

        input ->    - list of index pairs which are to be complemented together
                    - origin sequence string from GenBankParser

        output > complement string of nucleotides (TAGC)
        '''

        # sets checks for < or > to return reverse complement instead of normal complement to False
        # as standard
        # reverse_complement = True

        nucleotide_dictionary = {'A':'t', 'T':'a', 'C':'g', 'G':'c'}

        # cleans index pair string from 'complement' and the brackets '(', ')'
        cleaned_index_pair_string = index_pair_string.replace('complement(', '').replace(')', '')

        # splits the cleaned index string on the double dot (..)
        split_index_pair_list = cleaned_index_pair_string.split('..')

        # gets the start index from the split index pairs
        start_index = split_index_pair_list[0]

        # gets the end index from the split index pairs
        end_index = split_index_pair_list[1]

        # print('start is ', start_index, 'end is ', end_index)
        #The following block checks the start and stop indices for markers (> or <)
        # checks for a '>' in the start index
        if '>' in start_index:
            # reverse_complement = True
            start_index = int(start_index.replace('>', ''))

        # checks for a '<' in the start index
        elif '<' in start_index:
            # reverse_complement = True
            start_index = int(start_index.replace('<', ''))

        else:
            start_index = int(start_index)

        # checks for '>' in the end index
        if '>' in end_index:
            # reverse_complement = True
            end_index = int(end_index.replace('>', ''))

        # checks for the '<' in end index
        elif '<' in end_index:
            # reverse_complement = True
            end_index = int(end_index.replace('<', ''))

        else:
            end_index = int(end_index)

        if self.mode_select == False:
            # this list comprehension creates a sequence specified by the indices from the features
            # file and joins it together
            nucleotide_string = GenBankParser.origin_sequence[start_index-1:end_index]

            # returns a string containning the complement of the input string.
            # loops over the input string and joins the value of each key from the dictionary.
            # i is every seperate nucleotide in the sequence. .upper() makes sure it is in uppercase
            complement_str = ''.join([nucleotide_dictionary[i] for i in nucleotide_string.upper()])

            # checks if the reverse complement marker is True and inverts the sequence if so.
            # if reverse_complement == True:
            #     return complement_string[::-1]

            return complement_str[::-1]

        elif self.mode_select == True:
            upper_nuc_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
            # this list comprehension creates the sequence specified by the indices from a features
            # file and joins it together
            temp_string = GenBankParser.origin_sequence[:start_index-1]
            # this list comprehension creates a sequence specified by the indices from the features
            # file and joins it together, takes the nucleotides from the lowercase dict
            complement_string = ''.join([nucleotide_dictionary[i] for i in temp_string.upper()])

            # reassigns the temp_string with the wanted sequence indices
            temp_string = GenBankParser.origin_sequence[start_index-1:end_index].upper()

            # returns a string containning the complement of the input string.
            # loops over the input string and joins the value of each key from the dictionary.
            # i is every seperate nucleotide in the sequence. .upper() makes sure it is in uppercase
            complement_string += ''.join([upper_nuc_dict[i] for i in temp_string.upper()])
            # complement_string += temp_string
            # Returns the reverse complement
            return complement_string[::-1]


    def join_string(self, index_pair_string):
        '''
        This function loops over a given list of indices from the sequence and joins the sequence
        chuncks together. Works for uppercase as well.

        if uppercase is True

        In this statement the mode is in uppercase so the join works a little different.
        The sequence between each joined sequence chunck is in lowercase while the actual
        joined sequences are in uppercase. the piece of code in this part of the elif will
        generate the right output by splitting the sequence index pairs that are next to each
        -other.

        input ->    - list of index pairs which are to be joined together
                    - origin sequence string from GenBankParser

        ouput ->    - joined sequence junks as a string
        '''
        # take every index pair from the list and split them on the double dot (..)
        output_seq = ''

        # remove the .join, the opening bracket: (, and the closing bracket: )
        cleaned_index_string = index_pair_string.replace('join(', '').replace(')', '')

        # Split the seperate index pairs on the comma (,)
        split_cleaned_index_list = cleaned_index_string.split(',')


        # check the mode selector for seperated mode
        if self.mode_select == False:

            # loops over the list of index pairs (start_index..end_index)
            for index_pair in split_cleaned_index_list:
                #   check if index in pairs or single index
                if '..' in index_pair:
                #   seperate start_index from end_index
                    seperated_indices = index_pair.split('..')

                #   takes the first digit (start) and converts it to an integer
                    start_index = int(seperated_indices[0])

                #   takes the second integer (end) and converts it to an integer
                    end_index = int(seperated_indices[1])

                    # adds each sequence chunck to the output string i.e. 'joins' them
                    output_seq += GenBankParser.origin_sequence[start_index-1:end_index].strip()
                else:
                    start_index = int(index_pair)
                    output_seq += GenBankParser.origin_sequence[start_index-1].strip()

        # checks the mode selector for uppercase mode
        elif self.mode_select == True:

            # gets amount of index_pairs / indices
            length = len(split_cleaned_index_list)
            index_pair_index = 0

            start_index_sequence = split_cleaned_index_list[0].split('..')
            start_index_for_lowercase = int(start_index_sequence[0])-1

            # adds the sequence in lowercase from the start of the sequence until the start index-1
            # for the first sequence parts to be joined.
            output_seq = GenBankParser.origin_sequence[:start_index_for_lowercase]

            while index_pair_index < length-1: # loops over the list until the last entry

                # loops over line index and increments by 1 each loop
                for line_index in range(index_pair_index,index_pair_index+1):

                    # check if the first index pair is a pair or a single number
                    if '..' in split_cleaned_index_list[line_index]:
                        # get the first index pair and split it on the (..)
                        first_index_pair = split_cleaned_index_list[line_index].split('..')
                        # get the start index from the first pair
                        start_first = int(first_index_pair[0])-1
                        # get the end index from the first pair
                        end_first = int(first_index_pair[1])
                        # adds the required parts of join to the output sequence in uppercase
                        output_seq += GenBankParser.origin_sequence[start_first:end_first].upper()
                    else:
                        start_first = int(first_index_pair[0])
                        output_seq += GenBankParser.origin_sequence[start_first].upper()

                    # checks if second index is a pair or a singlenumber
                    if '..' in split_cleaned_index_list[line_index+1]:
                        # split the second index pair
                        second_index_pair = split_cleaned_index_list[line_index+1].split('..')
                        # second block ends one index early so no double letters in sequence.
                        start_second = int(second_index_pair[0])-1

                    else:
                        # second block ends one index early so no double letters in sequence.
                        start_second = int(split_cleaned_index_list[line_index+1])-1

                    # adds the parts between the joined sequence to the next joined sequence in
                    # lower case
                    output_seq += GenBankParser.origin_sequence[end_first:start_second]

                    # increments the line_index by 1
                    index_pair_index += 1

            # adds in the last join manually because it will be skipped otherwise

            last_index_pair = split_cleaned_index_list[index_pair_index]
            if '..' in last_index_pair:
                last_index_pair = split_cleaned_index_list[index_pair_index].split('..')
                start_index = int(last_index_pair[0])
                end_index = int(last_index_pair[1])
                output_seq += GenBankParser.origin_sequence[start_index-1:end_index].upper()
            else:
                start_index = int(last_index_pair)
                output_seq += GenBankParser.origin_sequence[start_index].upper()


        return output_seq


    def return_sequence(self, index_string):
        '''
        This function takes the start and stop index of the feature and returns the corresponding
        chunck of the sequence. This is the standard index to sequence converted, i.e. without any
        join or complement statements.

        input ->    - index_string = string containing the start and end indices as (1..2)

        output ->   - sequence in the form of a string
        '''

        # checks if the index string contains one or two indices and handles accordingly
        if '..' in index_string:
            # splits the two index numbers (strart:stop) and puts them into a temp_list
            temp_index_list = index_string.split('..')
            # takes the first digit (start) and converts it to an integer
            start_index = int(temp_index_list[0])-1
            # takes the second integer (end) and converts it to an integer
            end_index = int(temp_index_list[1])
            if self.mode_select == False:
                return GenBankParser.origin_sequence[start_index:end_index]
            elif self.mode_select == True:
                temp_string = GenBankParser.origin_sequence[:start_index]
                temp_string += GenBankParser.origin_sequence[start_index:end_index].upper()
                return temp_string
        else:
            start_index = int(index_string)-1
            if self.mode_select == False:
                return GenBankParser.origin_sequence[start_index]
            elif self.mode_select == True:
                temp_string = GenBankParser.origin_sequence[:start_index]
                temp_string += GenBankParser.origin_sequence[start_index].upper()
                return temp_string


    def complement_join(self, sequence_string):
        '''
        This function runs when a complement(join(index_numbers)) is encountered.
        It takes the sequence put out by the join function and complements it.

        input ->    - string of joined sequences

        output ->   - complemented string (atcg -> tagc)

        '''
        if self.mode_select == False: # if selected mode is 'seperated'
            nucleotide_dictionary = {'A':'t', 'T':'a', 'C':'g', 'G':'c'}
            complement_string = ''.join([nucleotide_dictionary[i] for i in sequence_string.upper()])
        elif self.mode_select == True:
            # this dictionairy can distinguish between upper and lowercase input.
            upper_nuc_dict = {'a':'t','A':'T','t':'a','T':'A','c':'g','C':'G','g':'c','G':'C'}
            # returns a string containning the complement of the input string.
            # loops over the input string and joins the value of each key from the dictionary.
            # i is every seperate nucleotide in the sequence. .upper() makes sure it is in uppercase
            complement_string = ''.join([upper_nuc_dict[i] for i in sequence_string])
        # return reverse of the complement string
        return complement_string[::-1]


    def order_amino_acids(self, sequence_string):
        '''
        This functions takes a string of amino acids and returns it sorted alphabetically.

        input:  - string of amino acids

        output: - sorted/ordered string of amino acids (a to z)
        '''
        sequence_string = sequence_string
        return ''.join(sorted(sequence_string))


    def sequence_source_isolator(self):
        '''
        This function isolates the /gene= or /source= parts of the feature and returns it as string

        input ->    -self.features > list containing features line

        output ->   -isolated source string
        '''
        # Take the first entry of the list [1] ('gene="CFTR"\n') and return it.
        return str(self.features[1]).strip()


    def name_isolator(self):
        """
        This function isolates the name of the feature and returns the name as a string.

        input ->    - self.features = a list containing the lines from the file as a list

        output ->   - name as a string.
        """
        # Take the first entry of the list [0] (example:'exon 87858..88040\n '), split it in parts.
        split_first_line = self.features[0].split()
        # from the list with the split entries, take the first entry, which is the name.
        name_container = split_first_line[0]
        # print(name_container)
        return name_container


    def string_length_modifier(self, sequence_string):
        '''
        This function slices a continuous sequence string up into lines of 60 characters and adds a
        new line(\\n). Then adds the following line on the next line.

        - parameter : Sequence string

        -> multiple lines with a max length of 60 characters.
        '''
        new_string = ''
        max_line_length = 60
        for sequence_index in range(0, len(sequence_string), max_line_length):
            new_string += sequence_string[sequence_index:sequence_index+max_line_length] + '\n'
        return new_string


    def finished_output(self):
        '''
        This function puts all the different parts together and formats it correctly for the
        output file
        '''
        try:
            # makes a string in the format: '>' <feature_name> <feature_source>
            output_string = '>' + self.features[0].split()[0] +' /'+ self.sequence_source.strip() +'\n'
            output_string += self.output_sequence
            return output_string
        except TypeError:
            print('TypeError while finishing output')




def main():
    '''
    standard main function
    '''
    Features

if __name__ == '__main__':
    main()
