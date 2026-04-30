import argparse
import os
import socket

from huggingface_hub import HfApi, create_repo


def force_ipv4_if_requested():
    if os.environ.get("HF_FORCE_IPV4") != "1":
        return
    original_getaddrinfo = socket.getaddrinfo

    def getaddrinfo_ipv4(*args, **kwargs):
        return [
            item for item in original_getaddrinfo(*args, **kwargs)
            if item[0] == socket.AF_INET
        ]

    socket.getaddrinfo = getaddrinfo_ipv4


def parse_args():
    parser = argparse.ArgumentParser(description="Upload a DeepARG Transformers bundle to Hugging Face Hub.")
    parser.add_argument("--repo-id", required=True, help="Hub repo id, for example username/deeparg-ls.")
    parser.add_argument("--folder", required=True, help="Local folder produced by save_pretrained/export.")
    parser.add_argument("--private", action="store_true")
    parser.add_argument("--repo-type", default="model")
    parser.add_argument("--commit-message", default="Upload DeepARG Transformers model")
    parser.add_argument("--token", default=None, help="Optional HF token; otherwise uses cached login/HF_TOKEN.")
    return parser.parse_args()


def main():
    force_ipv4_if_requested()
    args = parse_args()
    create_repo(
        repo_id=args.repo_id,
        repo_type=args.repo_type,
        private=args.private,
        token=args.token,
        exist_ok=True,
    )
    api = HfApi(token=args.token)
    api.upload_folder(
        repo_id=args.repo_id,
        repo_type=args.repo_type,
        folder_path=args.folder,
        commit_message=args.commit_message,
    )


if __name__ == "__main__":
    main()
