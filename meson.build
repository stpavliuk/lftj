project('wcoj', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++17'])

lftj_lib = static_library(
  'lftj_core',
  ['trie.cpp'],
)

wcoj = executable(
  'wcoj',
  'main.cpp',
  link_with : lftj_lib,
  install : true,
)

lftj_tests = executable(
  'lftj_tests',
  'test_lftj.cpp',
  link_with : lftj_lib,
)

test('lftj-tests', lftj_tests)
