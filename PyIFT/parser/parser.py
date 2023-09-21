#!/usr/bin/env python3

import os
import re
import argparse
from os.path import join as pjoin
from os.path import basename


def add_pointer_function(out_file, stream, func_dict):
    global args
    if use_stable(args) and 'stable' not in func_dict.keys():
        return

    func_header = stream.split(';')[0]
    # regex to remove everything except function name
    func_name = re.sub(r'(^\s*[a-zA-Z0-9\_\*\s]*\s\**)|(\([\,\.a-zA-Z0-9\s\*\_]*\))', '', func_header)

    exit_if_stable_error(args, func_name)
    try_write_summary(func_name)

    buffer = '%newobject ' + func_name + ';\n' + autodoc + func_header + ';\n\n'
    out_file.write(buffer)
    global func_count
    func_count += 1


def add_pointer_pyfunc(out_file, stream, func_dict):
    global args
    if use_stable(args) and 'stable' not in func_dict.keys():
        return

    func_header = stream.split('{')[0]
    # regex to remove everything except function name
    func_name = re.sub(r'(^\s*[a-zA-Z0-9\_\*]*\s\**)|(\([\,\.a-zA-Z0-9\s\*\_]*\))\s*', '', func_header)

    exit_if_stable_error(args, func_name)
    try_write_summary(func_name)

    buffer = '%newobject ' + func_name + ';\n' + autodoc + func_header + ';\n\n'
    out_file.write(buffer)
    global func_count
    func_count += 1


def add_normal_function(out_file, stream, func_dict):
    global args
    if use_stable(args) and 'stable' not in func_dict.keys():
        return

    func_header = stream.split(';')[0]
    if '{' in func_header:
        # when function is declared on header
        func_header = stream.split('{')[0]

    func_name = re.sub(r'(^\s*[a-zA-Z0-9\_\*\s]*\s\**)|(\([\,\.a-zA-Z0-9\s\*\_]*\))\s*', '', func_header)

    exit_if_stable_error(args, func_name)
    try_write_summary(func_name)

    buffer = autodoc + func_header + ';\n\n'
    out_file.write(buffer)
    global func_count
    func_count += 1


def add_normal_pyfunc(out_file, stream, func_dict):
    global args
    if use_stable(args) and 'stable' not in func_dict.keys():
        return

    func_header = stream.split('{')[0]
    func_name = re.sub(r'(^\s*[a-zA-Z0-9\_\*]*\s\**)|(\([\,\.a-zA-Z0-9\s\*\_]*\))\s*', '', func_header)

    exit_if_stable_error(args, func_name)
    try_write_summary(func_name)

    buffer = func_header + ';\n\n'
    out_file.write(buffer)
    global func_count
    func_count += 1


def add_destroyer(out_file, object, function):
    buffer = '\t~' + object + '() {\n' \
             '\t\t' + object + '* ptr = ($self);\n' \
             '\t\t' + function + '(&ptr);\n' \
             '\t}\n'
    out_file.write(buffer)


def add_write(out_file, object, function):
    buffer = '\tvoid Write(char* filename) {\n' \
             '\t\t' + object + '* ptr = ($self);\n' \
             '\t\t' + function + '(ptr, filename);\n' \
             '\t}\n'
    out_file.write(buffer)


def add_extension(out_file, ext_file):
    ext_path = pjoin(ift_dir, 'PyIFT', 'extend', ext_file)
    with open(ext_path, 'r') as file:
        for line in file.read().split('\n'):
            out_file.write('\t' + line + '\n')
        out_file.write('\n')


def parse_struct_extend(out_file, stream, func_dict):

    typedef = re.findall(r'^\s*[a-zA-Z0-9]+\s+[a-zA-Z0-9\_]+', stream)[0]
    inside = re.findall(r'\{[a-zA-Z0-9\s\;\:\,\_\/\+\-\=\\\*\@\]\[\)\(\.\>\<\'\&]+\}[a-zA-Z0-9\s\_]*\;', stream)[0]

    if 'name' in func_dict.keys():
        # when struct is typedefined on other part of the file
        object = func_dict['name']
        out_file.write('\n\ntypedef struct ' + inside[:-1] + ' ' + object + ';\n\n')
    else:
        # when struct is created with typedef
        object = re.findall(r'\}\s*[a-zA-Z0-9]+\s*(?=\;)', inside)[0]
        object = re.sub(r'\}\s*', '', object)
        out_file.write('\n\n' + typedef + ' ' + inside + '\n\n')

    # begin swig extend
    out_file.write('%extend ' + object + ' {\n\n')

    for key in func_dict.keys():
        if key == 'destroyer':
            add_destroyer(out_file, object, func_dict[key])
        elif key == 'extend':
            add_extension(out_file, func_dict[key])
        # elif key == 'write':
        #   add_write(out_file, object, func_dict[key])

    # end extend
    out_file.write('};\n\n')
    global object_count
    object_count += 1


def parse_enum(out_file, stream, func_dict):

    typedef = re.findall(r'^\s*[a-zA-Z0-9]+\s+[a-zA-Z0-9]+', stream)[0]
    inside = re.findall(r'\{[a-zA-Z0-9\s\;\:\,\_\/\+\-\\\*\@\]\[\)\(\.\=\|]+\}[a-zA-Z0-9\s\_]+\;', stream)[0]

    out_file.write('\n\n' + typedef + ' ' + inside + '\n\n')


def find_includes(content):
    updated = re.sub(r'\#include', '%include', re.sub(r'\.h', '.i', content))
    out = re.findall(r'\%include\ *\"ift[a-zA-Z]+\.i\"', updated)
    return out


def find_arguments(chunk, file_name):
    out = {}
    arguments = re.findall(r'[a-zA-Z0-9\s\=\,\.]*(?=\))', chunk)
    if arguments != []:
        for arg in re.split(r'\,\s*', arguments[0]):
            tmp = re.split(r'\s*\=\s*', arg)
            key = tmp[0]
            if len(tmp) > 1:
                out[key] = tmp[1]
            else:
                out[key] = ''
            if key != '' and key not in argument_list:
                print('\nATTENTION! ', key, ' is not a swig export argument, file', file_name)
    return out


def exit_if_stable_error(args, func_name):
    if use_stable(args) and is_stable_ok(func_name) is False:
        print('Stable function %s is not available on the demo directory' % func_name)
        print('PyIFT parsing stopped!')
        exit(-1)


def is_stable_ok(func):
    return True # when Falcao requested the additional functions to be release for the image analysis course
                # I removed the necessity of functions to be used on the demo directory
    # global demo_glob
    # name = re.sub(r'ift', '', func)
    # return bool(re.search(name, demo_glob))


def use_stable(args):
    return args.stable is not None


def build_summary(args):
    return args.summary is not None


def try_write_summary(func_name):
    global sum_file
    if sum_file:
        sum_file.write("\t\t" + func_name + "\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="PyIFT parser")
    parser.add_argument('--stable', help='The PyIFT parser will select only the stable funcions which '
                                           'contains a demo file', required=False)
    parser.add_argument('--summary', help='Create summary.txt file showing which functions from which files'
                                          'were added by the parser', required=False)

    args, _ = parser.parse_known_args()

    func_count = 0
    object_count = 0
    autodoc = '%feature("autodoc", "2");\n'
    ift_dir = os.environ['NEWIFT_DIR']
    demo_glob = ''

    sum_file = None
    if build_summary(args):
        sum_file = open("summary.txt", "w")

    if use_stable(args):
        for dirpath, dirnames, filenames in os.walk(pjoin(ift_dir, 'PyIFT', 'demo')):
            for filename in filenames:
                if '.py' in filename or '.ipynb' in filename:
                    with open(pjoin(dirpath, filename)) as f:
                        demo_glob = demo_glob + f.read()

    headers_dir = pjoin(ift_dir, 'include')
    headers_files = []
    for dirpath, dirnames, filenames in os.walk(headers_dir):
        for filename in [f for f in filenames if f.endswith(".h")]:
            headers_files.append(pjoin(dirpath, filename))

    swig_dir = pjoin(ift_dir, 'PyIFT', 'src')
    if not os.path.exists(swig_dir):
        os.makedirs(swig_dir)
    swig_file_list = []

    pyfunc_dir = pjoin(ift_dir, 'PyIFT', 'pyfunc')
    pyfunc_files = os.listdir(pyfunc_dir)

    argument_list = ['newobject', 'extend', 'destroyer', 'name', 'stable']

    if sum_file:
        sum_file.write("libIFT files\n")

    for header in headers_files:
        # because .cuh (cuda headers) have same name as C files
        if '.h' in header:

            file_name = basename(header).split('.')[0]
            swig_filepath = pjoin(swig_dir, (file_name + '.i'))
            swig_file_list.append(swig_filepath)

            if sum_file:
                sum_file.write("\n\t" + file_name + "\n")

            with open(pjoin(headers_dir, header), 'r', encoding="utf-8") as header_file:
                with open(swig_filepath, 'w') as swig_file:

                    file_content = header_file.read()

                    includes = find_includes(file_content)
                    for include in includes:
                        swig_file.write(include + '\n')
                    swig_file.write('\n')

                    maybe_export = re.split(r'\/\/\!\s*', file_content)

                    for chunk in maybe_export:

                        if 'swig(' == chunk[:5]:

                            export_args = find_arguments(chunk, file_name)

                            stream = chunk.split('\n', 1)[1]
                            next_line = chunk.split('\n', 2)[1] # deixar mais inteligente

                            if next_line == '':
                                print('Linebreak between swig notation and C code on file ', file_name,
                                      ' this might cause problems generating the swig files.')

                            if re.findall(r'(^|\s)struct(\s|\{)', next_line):
                                parse_struct_extend(swig_file, stream, export_args)
                            elif re.findall(r'(^|\s)enum(\s|\{)', next_line):
                                parse_enum(swig_file, stream, export_args)
                            else:
                                if 'newobject' in export_args.keys():
                                    add_pointer_function(swig_file, stream, export_args)
                                else:
                                    add_normal_function(swig_file, stream, export_args)


            with open(pjoin(ift_dir, 'PyIFT', 'pyift', 'pyift.i'), 'w') as pyift:

                pyift.write('%include "extend/iftSwig.i"\n\n')

                for include in swig_file_list:
                    pyift.write('%include "src/' + basename(include) + '"\n\n')

    if sum_file:
        sum_file.write("\n\nPyIFT exclusive files\n")

    for pyfunc in pyfunc_files:

        file_name = basename(pyfunc).split('.')[0]
        if sum_file:
            sum_file.write("\n\t" + file_name + "\n")

        with open(pjoin(pyfunc_dir, pyfunc), 'r') as pyfunc_file:
            with open(pjoin(swig_dir, pyfunc), 'a') as swig_file:

                file_content = pyfunc_file.read()

                swig_file.write('%inline %{\n\n')
                for line in file_content.split('\n'):
                    if line[:2] != '//':
                        swig_file.write(line + '\n')
                swig_file.write('\n')
                swig_file.write('%}\n\n')

                maybe_export = re.split(r'\/\/\!\s*', file_content)

                for chunk in maybe_export:

                    if 'swig(' == chunk[:5]:

                        export_args = find_arguments(chunk, pyfunc_file)

                        stream = chunk.split('\n', 1)[1]
                        next_line = chunk.split('\n', 2)[1] # deixar mais inteligente

                        if next_line == '':
                            print('Linebreak between swig notation and C code on file ', pyfunc_file,
                                  ' this might cause problems generating the swig files.')

                        if 'newobject' in export_args.keys():
                            add_pointer_pyfunc(swig_file, stream, export_args)
                        else:
                            add_normal_pyfunc(swig_file, stream, export_args)

    if sum_file:
        sum_file.close()

    print('\nWrapped ' + str(object_count) + ' objects and ' + str(func_count) + ' functions.\n')

