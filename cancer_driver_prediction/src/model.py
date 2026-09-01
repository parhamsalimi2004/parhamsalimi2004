import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import xgboost as xgb

def make_fake_feature_table(n_mutations: int=2000, seed=0):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({
        "damage_score": rng.uniform(0,1,n_mutations),
        "network_centrality": rng.uniform(0,1,n_mutations),
        "cooccurrence_score": rng.uniform(-1,1,n_mutations),
        "is_hotspot": rng.integers(0,2,n_mutations), #1=hotspot, 0=nonhotspot
        "label": rng.integers(0,2,n_mutations) #1=driver, 0=passenger
    })
    return df

def train_model(feature_table, feature_columns, label_column="label"):
    x = feature_table[feature_columns]
    y = feature_table[label_column]
    x_train, x_test, y_train, y_test = train_test_split(
        x, y, test_size=0.2, random_state=0, stratify=y)

    model = xgb.XGBClassifier(
        n_estimators=200, max_depth=4, learning_rate=0.05, eval_metric="logloss")
    model.fit(x_train,y_train)

    return model, x_test, y_test

def evaluate_model(model, x_test, y_test, hotspot_column="is_hotspot"):
    predictions = model.predict_proba(x_test)[:,1]
    overall_auc = roc_auc_score(y_test, predictions)
    print(f"overall AUC: {overall_auc:.3f}")

    hotspot_mask = x_test["is_hotspot"] == 1
    non_hotspot_mask = x_test["is_hotspot"] == 0
    hotspot_auc = roc_auc_score(y_test[hotspot_mask],predictions[hotspot_mask])
    non_hotspot_auc = roc_auc_score(y_test[non_hotspot_mask], predictions[non_hotspot_mask])

    print(f"hotspot AUC: {hotspot_auc: .3f}")
    print(f"non_hotspot AUC: {non_hotspot_auc:.3f}")

    return overall_auc

if __name__ == "__main__":
    print("Testing the pipeline plumbing with fake data (not real results)\n")

    fake_data = make_fake_feature_table()
    feature_cols = ["damage_score", "network_centrality", "cooccurrence_score", "is_hotspot"]

    model, X_test, y_test = train_model(fake_data, feature_cols)
    evaluate_model(model, X_test, y_test)
