import unittest

from rsmviewer.utils.logger import PluginLogger


class LoggerTests(unittest.TestCase):
    def test_debug_is_silent_by_default_and_enabled_explicitly(self):
        logger = PluginLogger()
        self.assertFalse(logger.debug_enabled)
        logger.set_debug(True)
        self.assertTrue(logger.debug_enabled)
        logger.set_debug(False)
        self.assertFalse(logger.debug_enabled)


if __name__ == "__main__":
    unittest.main()