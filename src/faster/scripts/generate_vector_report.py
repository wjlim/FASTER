#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
from ..reports.vector_report import VectorReport

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Generate STR vector analysis report')
    parser.add_argument('--vector-dir', required=True,
                      help='Directory containing vector.json files')
    parser.add_argument('--output-dir', required=True,
                      help='Directory to save report outputs')
    
    args = parser.parse_args()
    
    try:
        # Generate report
        report = VectorReport(args.vector_dir, args.output_dir)
        report_path = report.generate_report()
        logger.info(f"Generated vector report at: {report_path}")
        
    except Exception as e:
        logger.error(f"Error generating report: {str(e)}")
        raise

if __name__ == '__main__':
    main() 