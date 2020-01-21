#!/usr/bin/env python

"""Tests for `solexs_pipeline` package."""


import unittest
from click.testing import CliRunner

from solexs_pipeline import solexs_pipeline
from solexs_pipeline import cli


class TestSolexs_pipeline(unittest.TestCase):
    """Tests for `solexs_pipeline` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        assert 'solexs_pipeline.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
