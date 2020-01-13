pipeline {
  agent { 
    docker { 
      image 'conda/miniconda3-centos7:latest'
    } 
  }
  stages {
    stage('root-setup') {
      steps {
        slackSend (message: "Build Started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        sh 'yum install -y gcc gcc-c++ git gcc-gfortran libboost-dev make wget'
        dir('idaes-dev') {
          git url: 'https://github.com/makaylas/idaes-dev.git',
          credentialsId: '6ca01274-150a-4dd4-96ec-f0d117b0ea95'
        }
      }
    }
    stage('idaes-dev 3.7-setup') {
      steps {
        sh '''
         cd idaes-dev
         conda create -n idaes3.7 python=3.7 pytest
         source activate idaes3.7
         pip install -r requirements-dev.txt --user jenkins
         export TEMP_LC_ALL=$LC_ALL
         export TEMP_LANG=$LANG
         export LC_ALL=en_US.utf-8
         export LANG=en_US.utf-8
         python setup.py develop
         export LC_ALL=$TEMP_LC_ALL
         export LANG=$TEMP_LANG
         source deactivate
         '''
      }
    }
    stage('idaes-ext build') {
      steps {
        sh '''
         cd ..
         ls
         source activate idaes3.7
         bash scripts/compile_solvers.sh
         bash scripts/compile_libs.sh
         ls
         ls dist-lib
         ls dist-solvers
         idaes get-extensions
         source deactivate
         '''
      }
    }
    stage('idaes-dev 3.7-test') {
      steps {
        catchError(buildResult: 'SUCCESS', stageResult: 'FAILURE') {
          sh '''
           cd idaes-dev
           source activate idaes3.7
           pylint -E --ignore-patterns="test_.*" idaes || true
           pytest -c pytest.ini idaes -m "not nocircleci"
           source deactivate
           '''
        }
      }   
    }
  }
  post {
    success {
      slackSend (color: '#00FF00', message: "SUCCESSFUL - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    }

    failure {
      slackSend (color: '#FF0000', message: "FAILED - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
    }
  }
}
