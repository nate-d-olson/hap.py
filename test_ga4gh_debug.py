"""Test GA4GH module"""

import sys
sys.path.insert(0, 'src')

print("Starting ga4gh module import")

# Simple test - import one class directly  
try:
    print("Attempting to import GA4GHMetrics")
    from hap_py.haplo.ga4gh_compliance import GA4GHMetrics
    print("GA4GHMetrics imported successfully:", GA4GHMetrics)
except Exception as e:
    print("Failed to import GA4GHMetrics:", e)
    import traceback
    traceback.print_exc()

print("Module import complete")
