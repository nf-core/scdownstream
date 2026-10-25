"""Deterministic external-client replacement, loaded only by adapter.config."""


class CyteType:
    def __init__(self, adata, *, group_key, rank_key, auth_token):
        if rank_key not in adata.uns:
            raise ValueError("The adapter requires production rank data")
        self.adata = adata
        self.group_key = group_key

    def run(self, *, study_context, show_progress, timeout_seconds, require_artifacts):
        for field in ("annotation", "cellOntologyTerm", "cellOntologyTermID", "cellState"):
            self.adata.obs[f"cytetype_{field}_{self.group_key}"] = f"{self.group_key}:{field}"
        return self.adata
