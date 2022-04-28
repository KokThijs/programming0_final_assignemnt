'''
Name: genbankparser_class
This module contains my Features class to use in my main module.

Author: Thijs Kok
Student number: 348224
Class: Data Science for Life Science

This module contains the class that makes objects from each chunck of text.
'''
import re
import sys

class GenBankParser():
    '''
    The class GenBankParser is a class that creates an object by reading a given genbank file, it
    stores the file as a list and does operations on the list.
    It stores the 'origin', which is a sequence of DNA, mRNA or amino-acids.
    Furthermore it stores the trimmed definition given in the genbank file.
    '''
    origin_sequence = ''
    definition = ''
    isolated_feature = [] # list of isolated features
    list_of_feature_objects = []

    def __init__(self, input_file_name):
        # open the file
        try:
            with open(file=input_file_name, mode='r', encoding='utf-8') as input_file:
                self.file_as_list = input_file.readlines()
        except UnboundLocalError:
            sys.exit()
            # exception if the file does not exist
        except FileNotFoundError:
            print('File not found \n> please speficy a valid file in the cmd-line when running',
            'the script')
            sys.exit()

        # this block searches the line numbers of the search arguments in the input file
        origin_search_argument = 'ORIGIN'
        feature_search_argument = 'FEATURES'
        definition_search_argument = 'DEFINITION'
        self.origin_location = self.line_finder(origin_search_argument)
        self.features_location = self.line_finder(feature_search_argument)
        self.definition_location = self.line_finder(definition_search_argument)

        # assigns the list of indices from each feature to self
        self.features_index_list = self.features_lookup()

        # assigns the list of isolated features to self
        self.isolated_features = self.isolate_feats(self.features_index_list)

        # stores the origin sequence of the file as an attribute of all Genbankparser instances
        temp_string = self.sequence_isolater(self.file_as_list[self.origin_location::])
        GenBankParser.origin_sequence = temp_string

        # stores the definition of the file as an attribute of all Genbankparser instances
        GenBankParser.definition = self.definition_returner(self.file_as_list)

        # stores the list of isolated features of the file as an attribute of all Genbankparser
        # instances
        GenBankParser.isolated_feature = self.isolated_features


    def line_finder(self, search_argument):
        '''
        line_finder returns the line in which a given search argument is given.

        input:  - the content of the file, read into a list

                - string containing the argument that is looked for, for example ORIGIN

        output: index number of the line in which the search statement is found

        '''
        temp_list = []
        for line_number, line_in_list in enumerate(self.file_as_list, 0):
            # checksif the DEFINITION consists of two lines
            if search_argument in line_in_list and line_in_list.strip().endswith(','):
                temp_list.append(line_number)
                temp_list.append(line_number+1)
                return temp_list

            elif search_argument in line_in_list.rstrip():
                return line_number



    def definition_returner(self, file_as_list):
        '''
        This function returns the definition part of the input file as a string.

        input -> file as lines in list

        output -> Definition string
        '''
        if type(self.definition_location) == int:
            try:
                definition_string = file_as_list[self.definition_location] + '\n'
                # return stripped line-DEFINITION
                return definition_string.removeprefix('DEFINITION').strip()
            except TypeError:
                print('type_error')
        elif type(self.definition_location) == list:
            try:
                # start index = the first entry in the list, first line of the definition
                start_index = self.definition_location[0]
                # end index = the second entry in the list, second line of the definition
                end_index = self.definition_location[1]

                # store the string, remove 'DEFINITION' and remove any whitespace and \n
                temp_start = self.file_as_list[start_index].removeprefix('DEFINITION').strip() +' '
                # store the string and remove any whitespace and \n
                temp_end = self.file_as_list[end_index].strip()
                # concatenate both temporary strings
                definition_string = (temp_start + temp_end)
                return definition_string
            except IndexError:
                print('This index does not exist')


    def sequence_isolater(self, line_list):
        '''
        sequence_isolator isolates the the nucleotide sequence in the file under the header:'ORIGIN'
        in the input files the sequence is located from column 11 until 76. This function takes the
        characters that make up the sequence and puts them in one string.

        - input = list of every line in the read in file

        - output -> string containing the entire sequece
        '''
        self.seq_string = ''
        # reads individual lines from the list (starting from index one, so "origin" is skipped)
        for line in line_list[1::]:
            #reads individual character from each line from index 10 to 75
            for character in line:
                # looks for any entry containing a letter in lower- or uppercase
                if re.search(r"[a-zA-Z]", character):
                # adds each charachter to a list without whitespace and line endings
                    self.seq_string += character.strip()

        return self.seq_string


    def features_block_isolator(self, file_as_list):
        '''
        This function returns the block in the file containing all features.
        input:  - the list in which the file in read into

                - 'start_features' = index of the line in which the 'FEATURES' header is found.

                - 'end_features' = the index of the line in which the 'ORIGIN' header is found.

        output -> list containing only the features block
        '''
        # takes the lines from index start_features up to end_features (= start origin -1)
        try:
            isolated_features_block = file_as_list[self.features_location+1:self.origin_location]
        except IndexError:
            print('index out of range')
        return isolated_features_block


    def features_lookup(self): # needs list of features
        '''
        This Function looks for headers in the features part of the file and returns a list
        containing the respective line numbers.
        The function uses a regular expression that searches for a pattern containing 5 whitespaces
        and then a single letter.

        input:

        -   self.file_as_list > this is the entire file, read in as a list.

        output:

            - A list containing the indices where the regex pattern was found.
        '''
        pattern = r'^\s{5}[a-zA-Z]' # looks for 5 white_spaces and then a letter
        line_list = []
        line_number = 0 # = index from the list of features
        for line in self.file_as_list:
            try:
                if re.search(pattern, line): # checks if the regex pattern is in the current line

                # appends index at which the regex is found to a new list
                    line_list.append(line_number)
                else:
                    pass

                line_number += 1
            except:
                print('invalid feature')

        # print(line_list)
        return line_list


    def isolate_feats(self, feat_indices):
        '''
        This function takes the indices of where the feature headers are found with the features
        lookup function. It takes these these indices and makes a new list from one index to the
        next index and so on.

        input:  -list of indices of the features

        ouput:  -features in the form of individual lines in a list
        '''

        output_features_list = []
        length = len(feat_indices)# gets the amount of entries in the list of feature indices
        line_index = 0
        try:
            while line_index < length-1: # loops over the list until the last entry

                temp_data = ""
                for line in range(feat_indices[line_index], feat_indices[line_index+1]):
                    temp_data += self.file_as_list[line]

                output_features_list.append(''.join(temp_data))
                line_index += 1 # increments the line_index by 1

            # This adds the text from the last index in the list of feature indices until
            # the start of the origin block
            temp_data = ''
            for line in range(feat_indices[line_index],self.origin_location):
                temp_data += self.file_as_list[line]
            output_features_list.append(''.join(temp_data))
        except IndexError:
            print('index out of range 2')

        finally:
            return output_features_list


def main():
    '''
    standard main function, i'm not sure if it's nessecary but i just included it.
    '''
    GenBankParser


if __name__ == '__main__':
    main()
