pipeline {
  agent { 
    docker { 
      image 'conda/miniconda3-centos7:latest'
    } 
  }
  stages {
    // Commented out for an example until we start using these parameters
    // stage('cron-nightly-test') {
    //   when {
    //     expression { params.BUILD_SCHEDULE == 'Nightly'}
    //   }
    //   steps {
    //     sh 'echo "nightly works"'
    //   }
    // }
    // stage('cron-weekly-test') {
    //   when {
    //     expression { params.BUILD_SCHEDULE == 'Weekly'}
    //   }
    //   steps {
    //     sh 'echo "weekly works"'
    //   }
    // }

    // Until I get access to an HSL license, I'm just going to run the idaes-tests
    stage('root-setup') {
      steps {
        slackSend (message: "Build Started - ${env.JOB_NAME} ${env.BUILD_NUMBER} (<${env.BUILD_URL}|Open>)")
        sh 'yum install -y gcc g++ git gcc-gfortran libboost-dev make'
        sh 'cd ..'
        sh 'git clone git@github.com:makaylas/idaes-dev.git'
        sh 'cd idaes-dev'
      }
    }
    stage('3.7-setup') {
      steps {
        sh '''
         conda create -n idaes3.7 python=3.7 pytest
         source activate idaes3.7
         pip install -r requirements-dev.txt --user jenkins
         export TEMP_LC_ALL=$LC_ALL
         export TEMP_LANG=$LANG
         export LC_ALL=en_US.utf-8
         export LANG=en_US.utf-8
         python setup.py develop
         idaes get-extensions
         export LC_ALL=$TEMP_LC_ALL
         export LANG=$TEMP_LANG
         source deactivate
         '''
      }
    }
    stage('3.7-test') {
      steps {
        catchError(buildResult: 'SUCCESS', stageResult: 'FAILURE') {
          sh '''
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
