"""
Name: Final assessment programming0 2021

Author: Thijs Kok
Student number: 348224
Class: Data Science for Life Science

Description:
    - Module that imitates the "GenBank Feature Extractor" feature on the website:
    http://www.bioinformatics.org/sms2/genbank_feat.html

    - Takes a GenBank file containing DNA/mRNA/Proteïn sequence + features and parses out a
    detailed description of the features in the file and the corresponding location in the origin
    sequence.

Arguments:
    Pass the name of the Genbank file in the cmdline i.e.:
        "mRNA.gb"
"""

import sys
from features_class import Features
from genbankparser_class import GenBankParser


def mode_select(possible_choices):
    '''
    This function asks the user to select a mode while running the program.
    The seperated mode returns a file with each sequence part isolated from the others
    The uppercase mode returns a file with all previous entries from the sequence in underscore
    and the important part of the sequence in uppercase.

    -   input = user input > can select between 1 and 2 for seperated or uppercase respectively.

    -   output = True or False
    '''
    tried_choices = 0
    mode_choice = ""
    while mode_choice not in possible_choices:
        if tried_choices == 3:
            print('Is it really that hard to type the correct input?')
            mode_choice = input("Choose mode [%s]: " % ", ".join(possible_choices))
            tried_choices += 1
        elif tried_choices == 5:
            print('No valid mode chosen, returning file in seperated mode \n')
            return 'sep'
        else:
            mode_choice = input("Choose mode [%s]: " % ", ".join(possible_choices))
            tried_choices += 1
    return mode_choice


def create_feature_object(feature_list):
    """
    this function thakes the list of isolated features from self.isolated_features
    and makes an object out of each list contained in that list.
    """
    objects_list = []
    for feature_string in feature_list:
        feature_object = Features(feature_string)
        objects_list.append(feature_object)
    return objects_list


def file_writer(input_file_name, list_of_features):
    '''
    This function takes the input file's name, removes the file type extension and adds
    '_output.txt' on the end of the name. This will be the name of the output file.
    I.e. > CFTR_DNA.gb will become : CFTR_DNA_output.txt
    '''
    output_name_temp = input_file_name.removeprefix('.\\').split('.')
    output_name = output_name_temp[0] + '_features.txt'
    with open(output_name, mode='w', encoding='utf-8') as output_file_handler:
        output_file_handler.write(GenBankParser.definition + '\n\n')
        for feature_string in list_of_features:
            try:
                output_file_handler.write(feature_string.output_string +'\n')
            except TypeError:
                print('TypeError jammer dan 2')
    print(f'Feature extraction succesfull. Output written to:\n\n -> \t {output_name} \n')

def main():
    '''
    This is the main function that runs the module and classes.
    Here the mode select is ran and the file is returned accordlingly
    '''


    try:
        input_file_name = sys.argv[1]
        # checks the file for the right file extension (.gb or .gp)
        if input_file_name.endswith('.gb') or input_file_name.endswith('.gp'):
            try:
                # ask the user which mode he/she likes to use
                print('Please select a mode. sep for seperated mode, up for uppercase mode')
                mode_choice = mode_select(["sep", 'up'])
                if 'sep' in mode_choice:
                    print('\nSeperated mode selected, working on file... \n')
                    # leave the mode select at it's standard mode (False)
                    Features.mode_select = False
                elif 'up' in mode_choice:
                    print('\nUppercased mode selected, working on file... \n')
                    # select the uppercase mode (True)
                    Features.mode_select = True
            except:
                # if all fails, return the file in seperate mode
                print('\nsomething went wrong while selecting the mode',
                'returning file in seperated mode \n')
                Features.mode_select = False
            parser_object = GenBankParser(input_file_name)
            list_of_feature_objects = create_feature_object(parser_object.isolated_feature)
            GenBankParser.list_of_feature_objects = list_of_feature_objects
            # print(str(Features.object_count) +' '+ 'Features found an created')
            file_writer(input_file_name, list_of_feature_objects)
        else:
            print('Error, unexpected file extension.',
            'Please specify a valid genbank file with the \'.gb\' or \'.gp\' extension\n')
    except IndexError:
        print('Please speficy one genbank file (.gb or .gp) to open in the cmd-line')

if __name__ == '__main__':
    main()
