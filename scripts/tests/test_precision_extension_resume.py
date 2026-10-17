# Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
"""Pooled retained traces must not enter the uniform-budget resume writer."""
import asyncio
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'backend'))
from app.services import ground_validation as service


class PooledResumeTests(unittest.TestCase):
    def test_mixed_collection_is_rejected_before_source_or_trace_operations(self):
        job=dict(status='failed',precision_extension=dict(parent_id='preserved-study'),
                 request=dict(chains=24,n_measure=1000000))
        with patch.object(service,'get',return_value=job), \
             patch.object(service,'provenance',side_effect=AssertionError('No source operations allowed')), \
             patch.object(service,'save',side_effect=AssertionError('No retained-file writes allowed')):
            with self.assertRaisesRegex(RuntimeError,'frozen collector'):
                asyncio.run(service.resume('pooled-study'))

    def test_completed_original_study_keeps_its_existing_resume_rejection(self):
        with patch.object(service,'get',return_value=dict(status='complete')), \
             patch.object(service,'ACTIVE',{}):
            with self.assertRaisesRegex(RuntimeError,'Only stopped studies'):
                asyncio.run(service.resume('original-study'))


if __name__=='__main__':unittest.main()
