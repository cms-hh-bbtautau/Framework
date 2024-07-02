# FLAF

FLAF - Flexible LAW-based Analysis Framework.
Task workflow managed is done via [LAW](https://github.com/riga/law) (Luigi Analysis Framework).

## How to install
1. Setup ssh keys:
    - On GitHub [settings/keys](https://github.com/settings/keys)
    - On CERN GitLab [profile/keys](https://gitlab.cern.ch/-/profile/keys)

1. Clone the repository:
  ```sh
  git clone --recursive git@github.com:cms-flaf/Framework.git FLAF
  ```

1. Create a user customisation file `config/user_custom.yaml`. It should contain all user-specific modifications that you don't want to be committed to the central repository. Below is example of minimal content of the file (replace `USER_NAME` and `ANA_FOLDER` with your values):
    ```yaml
    fs_default:
        - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/'
    fs_anaCache:
        - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/'
    fs_anaTuple:
        - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/'
    fs_anaCacheTuple:
        - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/'
    fs_histograms:
    - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/histograms/'
    fs_json:
    - 'T3_CH_CERNBOX:/store/user/USER_NAME/ANA_FOLDER/jsonFiles/'
    analysis_config_area: config/HH_bbtautau
    compute_unc_variations: true
    store_noncentral: true
    ```

## How to load environment
1. Following command activates the framework environment:
    ```sh
    source env.sh
    ```

1. For the new installation or after you implement new law tasks, you need to update the law index:
    ```sh
    law index --verbose
    ```

1. Initialize voms proxy:
    ```sh
    voms-proxy-init -voms cms -rfc -valid 192:00
    ```

