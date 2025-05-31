node {
    stage 'Checkout'

    deleteDir()

    checkout scm

    stage 'Build / Install'
    sh '''
       cd ${WORKSPACE}
       pip install -e .
       pip install -r requirements-dev.txt
       pip install nox
    '''

    stage 'CI'
    sh '''
       cd ${WORKSPACE}
       nox -s lint format type_check tests
    '''
}
