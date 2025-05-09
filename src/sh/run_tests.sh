#!/bin/bash

##############################################################
# Test setup
##############################################################

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

# Initialize variables to track test failures
FAILED_TESTS=""
TEST_COUNT=0
FAIL_COUNT=0

##############################################################
# Boost unit tests
##############################################################

# echo "Running BOOST_TEST tests."

bin/test_haplotypes
TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Boost Unit tests - bin/test_haplotypes"
	FAILED_TESTS="${FAILED_TESTS}\n- Boost Unit tests (bin/test_haplotypes)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Boost Unit tests - bin/test_haplotypes"
fi

##############################################################
# Test Multimerge
##############################################################

/bin/bash ${DIR}/run_multimerge_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Multimerge tests - ${DIR}/run_multimerge_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Multimerge tests (${DIR}/run_multimerge_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Multimerge tests - ${DIR}/run_multimerge_test.sh"
fi

##############################################################
# Test Hapenum
##############################################################

/bin/bash ${DIR}/run_hapenum_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Hapenum test - ${DIR}/run_hapenum_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Hapenum test (${DIR}/run_hapenum_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Hapenum test - ${DIR}/run_hapenum_test.sh"
fi

##############################################################
# Test Hapcmp
##############################################################

/bin/bash ${DIR}/run_hapcmp_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Hapcmp test - ${DIR}/run_hapcmp_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Hapcmp test (${DIR}/run_hapcmp_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Hapcmp test - ${DIR}/run_hapcmp_test.sh"
fi

##############################################################
# Test Hap.py + path traversals
##############################################################

/bin/bash ${DIR}/run_pathtraversal_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Path traversal test - ${DIR}/run_pathtraversal_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Path traversal test (${DIR}/run_pathtraversal_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Path traversal test - ${DIR}/run_pathtraversal_test.sh"
fi

##############################################################
# Test Hap.py + FP regions
##############################################################

/bin/bash ${DIR}/run_fp_accuracy_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - FP region test - ${DIR}/run_fp_accuracy_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- FP region test (${DIR}/run_fp_accuracy_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - FP region test - ${DIR}/run_fp_accuracy_test.sh"
fi

##############################################################
# Test Hap.py + Faulty input variants
##############################################################

/bin/bash ${DIR}/run_faulty_variant_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Faulty variant test - ${DIR}/run_faulty_variant_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Faulty variant test (${DIR}/run_faulty_variant_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Faulty variant test - ${DIR}/run_faulty_variant_test.sh"
fi

##############################################################
# Test Hap.py safe leftshifting
##############################################################

/bin/bash ${DIR}/run_leftshift_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Leftshift test - ${DIR}/run_leftshift_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Leftshift test (${DIR}/run_leftshift_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Leftshift test - ${DIR}/run_leftshift_test.sh"
fi

##############################################################
# Test Hap.py + other VCF items
##############################################################

/bin/bash ${DIR}/run_other_vcf_tests.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Other VCF tests - ${DIR}/run_other_vcf_tests.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Other VCF tests (${DIR}/run_other_vcf_tests.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Other VCF tests - ${DIR}/run_other_vcf_tests.sh"
fi

##############################################################
# Test hom-ref block expansion and calls-only preprocessing
##############################################################

/bin/bash ${DIR}/run_gvcf_homref_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - GVCF hom-ref test - ${DIR}/run_gvcf_homref_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- GVCF hom-ref test (${DIR}/run_gvcf_homref_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - GVCF hom-ref test - ${DIR}/run_gvcf_homref_test.sh"
fi

##############################################################
# Test Hap.py + chr prefix detection
##############################################################

/bin/bash ${DIR}/run_chrprefix_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Chr prefix detection tests - ${DIR}/run_chrprefix_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Chr prefix detection tests (${DIR}/run_chrprefix_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Chr prefix detection tests - ${DIR}/run_chrprefix_test.sh"
fi

##############################################################
# Test Hap.py variant decomposition into primitives
##############################################################

/bin/bash ${DIR}/run_decomp_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Variant decomposition test - ${DIR}/run_decomp_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Variant decomposition test (${DIR}/run_decomp_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Variant decomposition test - ${DIR}/run_decomp_test.sh"
fi

##############################################################
# Test contig length calculation
##############################################################

${PYTHON} ${DIR}/run_fastasize_test.py
TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
    echo "[FAILED] - Contig length calculation test - ${DIR}/run_fastasize_test.py"
    FAILED_TESTS="${FAILED_TESTS}\n- Contig length calculation test (${DIR}/run_fastasize_test.py)"
    FAIL_COUNT=$((FAIL_COUNT + 1))
else
    echo "[PASSED] - Contig length calculation test - ${DIR}/run_fastasize_test.py"
fi

##############################################################
# Test Hap.py + integration
##############################################################

/bin/bash ${DIR}/run_integration_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Integration test - ${DIR}/run_integration_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Integration test (${DIR}/run_integration_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Integration test - ${DIR}/run_integration_test.sh"
fi

##############################################################
# Test Hap.py scmp allele and distance-based comparison
##############################################################

/bin/bash ${DIR}/run_scmp_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - SCMP test - ${DIR}/run_scmp_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- SCMP test (${DIR}/run_scmp_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - SCMP test - ${DIR}/run_scmp_test.sh"
fi

##############################################################
# Test Hap.py on tricky test cases
##############################################################

/bin/bash ${DIR}/run_giab_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Tricky indel test - ${DIR}/run_giab_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Tricky indel test (${DIR}/run_giab_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Tricky indel test - ${DIR}/run_giab_test.sh"
fi

##############################################################
# Test Performance + Consistency
##############################################################

/bin/bash ${DIR}/run_performance_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Performance / Consistency test - ${DIR}/run_performance_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Performance / Consistency test (${DIR}/run_performance_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Performance / Consistency test - ${DIR}/run_performance_test.sh"
fi

##############################################################
# Test GA4GH quantification
##############################################################

/bin/bash ${DIR}/run_quantify_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Quantify integration test - ${DIR}/run_quantify_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Quantify integration test (${DIR}/run_quantify_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Quantify integration test - ${DIR}/run_quantify_test.sh"
fi

##############################################################
# Test GA4GH stratified quantification
##############################################################

/bin/bash ${DIR}/run_quantify_stratification_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Quantify stratification test - ${DIR}/run_quantify_stratification_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Quantify stratification test (${DIR}/run_quantify_stratification_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Quantify stratification test - ${DIR}/run_quantify_stratification_test.sh"
fi


##############################################################
# Test PG Counting
##############################################################

/bin/bash ${DIR}/run_happy_pg_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - PG integration test - ${DIR}/run_happy_pg_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- PG integration test (${DIR}/run_happy_pg_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - PG integration test - ${DIR}/run_happy_pg_test.sh"
fi

##############################################################
# Test Hap.py + integration
##############################################################

/bin/bash ${DIR}/run_sompy_test.sh

TEST_COUNT=$((TEST_COUNT + 1))
if [[ $? -ne 0 ]]; then
	echo "[FAILED] - Som.py test - ${DIR}/run_sompy_test.sh"
	FAILED_TESTS="${FAILED_TESTS}\n- Som.py test (${DIR}/run_sompy_test.sh)"
	FAIL_COUNT=$((FAIL_COUNT + 1))
else
	echo "[PASSED] - Som.py test - ${DIR}/run_sompy_test.sh"
fi

##############################################################
# Test Summary
##############################################################
echo ""
echo "===== TEST SUMMARY ====="
echo "Tests Run: ${TEST_COUNT}"
echo "Tests Passed: $((TEST_COUNT - FAIL_COUNT))"
echo "Tests Failed: ${FAIL_COUNT}"

if [[ $FAIL_COUNT -gt 0 ]]; then
    echo -e "\nFailed Tests:${FAILED_TESTS}"
    # Return failure exit code but after running all tests
    exit 1
fi
