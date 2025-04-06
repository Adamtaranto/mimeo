#!/usr/bin/env python
"""
Main entry point for the mimeo command line tool.

This module provides a unified command-line interface for the mimeo package,
routing subcommands to the appropriate modules based on the first argument:

- mimeo x: Run cross-species comparison using run_interspecies
- mimeo self: Run self-alignment analysis using run_self
- mimeo map: Run mapping operations using run_map
- mimeo filter: Run filtering operations using run_filter

This approach allows for a more intuitive command structure while maintaining
compatibility with the individual command modules.
"""

import sys
from importlib import import_module


def main():
    """
    Parse command-line arguments and dispatch to the appropriate module.

    The first argument should be one of: x, self, map, filter
    All subsequent arguments are passed to the respective module's main function.
    """
    # Check if at least one argument was provided
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)

    # Get the subcommand from the first argument
    subcommand = sys.argv[1]

    # Define mapping of subcommands to modules
    commands = {
        'x': 'mimeo.run_interspecies',
        'self': 'mimeo.run_self',
        'map': 'mimeo.run_map',
        'filter': 'mimeo.run_filter',
    }

    # Check if the subcommand is valid
    if subcommand not in commands:
        print(f"Error: Unknown command '{subcommand}'")
        print_usage()
        sys.exit(1)

    # Remove the first argument (the subcommand) before passing to the module
    # This makes it look like the module was called directly with its own arguments
    sys.argv = [sys.argv[0]] + sys.argv[2:]

    # Import and run the appropriate module's main function
    try:
        module = import_module(commands[subcommand])
        module.main()
    except ImportError as e:
        print(f'Error importing module {commands[subcommand]}: {e}')
        sys.exit(1)
    except Exception as e:
        print(f"Error running command '{subcommand}': {e}")
        sys.exit(1)


def print_usage():
    """Display usage information for the mimeo command."""
    print("""
Usage: mimeo <command> [options]

Commands:
  x       Run cross-species comparison
  self    Run self-alignment analysis
  map     Run genomic mapping
  filter  Run filtering operations

For command-specific help:
  mimeo <command> --help
""")


if __name__ == '__main__':
    main()
