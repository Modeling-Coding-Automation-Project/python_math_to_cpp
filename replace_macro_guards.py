"""
replace_macro_guards.py

1. Searches all .hpp files under the current directory (excluding folders that
   start with '.' and the 'external_libraries' folder), then replaces every
   identifier of the form  __WORD__  with  WORD_  .

   Replacement rule (hpp):
     Leading  __  ->  (removed)
     Trailing __  ->  _
   Example:  __BASE_UTILITY_HPP__  ->  BASE_UTILITY_HPP_

2. Searches all .yml files under the '.github' folder and replaces every
   compiler-define flag of the form  -D__WORD__  with  -DWORD_  .

   Replacement rule (yml):
     -D__  ->  -D
     Trailing __  ->  _
   Example:  -D__BASE_UTILITY_USE_STD_COPY__  ->  -DBASE_UTILITY_USE_STD_COPY_
"""

import os
import re


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_hpp_files(root_dir: str) -> list[str]:
    """Return a list of .hpp file paths found under *root_dir*.

    Directories whose name starts with '.' or equals 'external_libraries'
    are skipped entirely (and not descended into).
    """
    hpp_files: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Modify dirnames in-place so os.walk does not descend into them.
        dirnames[:] = [
            d for d in dirnames
            if not d.startswith('.') and d != 'external_libraries'
        ]
        for filename in filenames:
            if filename.endswith('.hpp'):
                hpp_files.append(os.path.join(dirpath, filename))
    return hpp_files


# Matches identifiers that start with __ and end with __, e.g. __WORD__
# The lookbehind / lookahead make sure the __ is truly at an identifier
# boundary (not in the middle of a longer identifier).
_PATTERN = re.compile(r'(?<![A-Za-z0-9_])__([A-Za-z0-9_]+)__(?![A-Za-z0-9_])')

# Matches compiler-define flags of the form -D__WORD__, e.g. -D__BASE_UTILITY_USE_STD_COPY__
_D_PATTERN = re.compile(r'-D__([A-Za-z0-9_]+)__(?![A-Za-z0-9_])')


def find_yml_files(github_dir: str) -> list[str]:
    """Return a list of .yml file paths found under *github_dir* (.github)."""
    yml_files: list[str] = []
    if not os.path.isdir(github_dir):
        return yml_files
    for dirpath, _dirnames, filenames in os.walk(github_dir):
        for filename in filenames:
            if filename.endswith('.yml'):
                yml_files.append(os.path.join(dirpath, filename))
    return yml_files


def replace_macros(content: str) -> str:
    """Replace all ``__WORD__`` occurrences with ``WORD_``."""
    return _PATTERN.sub(r'\1_', content)


def replace_d_macros(content: str) -> str:
    """Replace all ``-D__WORD__`` occurrences with ``-DWORD_``."""
    return _D_PATTERN.sub(r'-D\1_', content)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _process_files(
    files: list[str],
    transform,
    changed: list[str],
    unchanged: list[str],
) -> None:
    """Apply *transform* to each file and record results."""
    for filepath in files:
        with open(filepath, 'r', encoding='utf-8') as fh:
            original = fh.read()

        modified = transform(original)

        if modified != original:
            with open(filepath, 'w', encoding='utf-8') as fh:
                fh.write(modified)
            changed.append(filepath)
            print(f"[MODIFIED]  {filepath}")
        else:
            unchanged.append(filepath)
            print(f"[unchanged] {filepath}")


def main() -> None:
    root_dir = os.path.dirname(os.path.abspath(__file__))

    # --- .hpp files ---
    hpp_files = find_hpp_files(root_dir)
    if not hpp_files:
        print("No .hpp files found.")
    else:
        print(f"Found {len(hpp_files)} .hpp file(s):\n")
        for f in hpp_files:
            print(f"  {f}")
        print()

    # --- .yml files under .github ---
    github_dir = os.path.join(root_dir, '.github')
    yml_files = find_yml_files(github_dir)
    if not yml_files:
        print("No .yml files found under .github.")
    else:
        print(f"Found {len(yml_files)} .yml file(s):\n")
        for f in yml_files:
            print(f"  {f}")
        print()

    if not hpp_files and not yml_files:
        return

    changed: list[str] = []
    unchanged: list[str] = []

    _process_files(hpp_files, replace_macros, changed, unchanged)
    _process_files(yml_files, replace_d_macros, changed, unchanged)

    print(
        f"\nDone. {len(changed)} file(s) modified, {len(unchanged)} file(s) unchanged.")


if __name__ == '__main__':
    main()
