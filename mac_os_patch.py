import subprocess
import platform

install_tool = 'install_name_tool'

ld_tool = 'otool'

platform_dir = platform.system()


def calculate_depth(lib_name):
    return lib_name.count('/') - 4


def change_id_command(lib_name, lib_path):

    return [install_tool, '-id', '@rpath/'+lib_name, lib_path]


def copy_dylib_command(lib_path):

    return ['cp']


def add_rpath_command(lib_name):

    dir_depth = calculate_depth(lib_name)
    dotdots = '../'*dir_depth

    return [install_tool, '-add_rpath', '@loader_path/'+dotdots+'_lib/Darwin', lib_name]


def change_lib_command(lib_name, referenced_lib, referenced_lib_path):

    return [install_tool, '-change', referenced_lib_path, '@rpath/'+referenced_lib, lib_name]


def get_dependencies(lib_name):

    output = subprocess.check_output(['otool', '-L', lib_name], shell=True).decode().split('\n')[1:]
    output = [item.strip() for item in output]
    output.remove('')
    output = [item.split()[0] for item in output]

    output = [item for item in output if 'libSystem' not in item]
    output = [item for item in output if '@rpath' not in item]

    return output


def get_filename(lib_path):
    return lib_path.split('/')[-1]


def get_filename_and_paths(lib_list):

    names = [get_filename(lib) for lib in lib_list]
    return dict(zip(names, lib_list))


def collect_dependency_tree(current_lib_list, dependency_set=set(), recursion_depth=0):

    if recursion_depth > 3:
        # Too many libraries!
        raise Exception('Too many dependencies to handle!')

    if len(current_lib_list) is 0:
        return dependency_set

    dependencies = set()
    lib_names = {}
    for lib in current_lib_list:

        lib_paths = get_dependencies(lib)
        lib_names.update(get_filename_and_paths(lib_paths))

        dependencies = dependencies.union(lib_names.keys())

    current_lib_names = get_filename_and_paths(dependency_set).keys()

    additional_libs = dependencies.difference(current_lib_names)

    new_lib_list = [lib_names[new_lib] for new_lib in additional_libs]

    '''
    print('Current:', current_lib_list)
    print('Have to check:', additional_libs)
    print('Already Have:', dependency_set)
    print(recursion_depth)
    print()
    '''

    dependency_set = dependency_set.union(new_lib_list)

    recursion_depth = recursion_depth + 1

    return collect_dependency_tree(new_lib_list, dependency_set, recursion_depth)


def patch_extension(lib_path):

    print(add_rpath_command(lib_path))
    subprocess.call(add_rpath_command(lib_path))

    depends_on = get_dependencies(lib_path)

    for dependency in depends_on:
        depend_file_name = get_filename(dependency)
        print(change_lib_command(lib_path, depend_file_name, dependency))
        subprocess.call(change_lib_command(lib_path, depend_file_name, dependency))


def patch_dylib(lib_path):

    # A little tricky since we need to copy the library itself.
    dest_dir = 'climt/_lib/'+platform_dir+'/'
    file_name = get_filename(lib_path)

    dest_path = dest_dir+file_name

    print(['cp', lib_path, dest_path])
    print(['chmod', '+w', dest_path])
    subprocess.call(['cp', lib_path, dest_path])
    subprocess.call(['chmod', '+w', dest_path])

    print(change_id_command(file_name, dest_path))
    subprocess.call(change_id_command(file_name, dest_path))

    depends_on = get_dependencies(dest_path)

    for dependency in depends_on:
        depend_file_name = get_filename(dependency)
        print(change_lib_command(dest_path, depend_file_name, dependency))
        subprocess.call(change_lib_command(dest_path, depend_file_name, dependency))

    print(['chmod', '-w', dest_path])
    subprocess.call(['chmod', '-w', dest_path])


def modify_library(library_path):

    extension = library_path.split('.')[-1]

    if extension == 'so':  # cython extension
        patch_extension(library_path)

    if extension == 'dylib':  # shared library
        patch_dylib(library_path)


lib_list = subprocess.check_output(['find', '.', '-name', "*.so"], shell=True).split()

lib_list = [element.decode() for element in lib_list]

libs_to_modify = list(collect_dependency_tree(lib_list))
libs_to_modify.extend(lib_list)

for lib in libs_to_modify:
    print(lib)
    modify_library(lib)
    print()
