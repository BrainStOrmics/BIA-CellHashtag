"""
HITL (Human-in-the-Loop) 交互。

提供 CLI 交互式提示，让用户在关键时刻介入决策。
"""


def hitl_prompt(message: str, options: dict) -> str:
    """
    向用户展示 HITL 提示并获取响应。

    Args:
        message: 提示消息。
        options: 可选操作 {key: description}。

    Returns:
        用户选择的 option key。
    """
    print("\n" + "=" * 60)
    print("🔍 需要人工确认 (Human-in-the-Loop)")
    print("=" * 60)
    print(message)
    print("\n可选操作:")
    for key, desc in options.items():
        print(f"  [{key}] {desc}")
    print("=" * 60)

    while True:
        choice = input("\n请选择 [{}]: ".format("/".join(options.keys()))).strip().lower()
        if choice in options:
            return choice
        print(f"无效选择，请输入: {', '.join(options.keys())}")


def parse_hitl_response(choice: str, context: dict) -> dict:
    """
    解析 HITL 响应，返回决策。

    Args:
        choice: 用户选择的 option key。
        context: 当前上下文。

    Returns:
        决策字典。
    """
    return {
        "action": choice,
        "context": context,
    }


def build_clustering_hitl_message(quality_results: dict, params: dict) -> str:
    """构建聚类质量评估的 HITL 消息。"""
    msg = "**聚类质量评估需要确认**\n\n"
    msg += f"当前参数:\n"
    msg += f"- Resolution: {params.get('resolution', 'N/A')}\n"
    msg += f"- 聚类数: {params.get('n_clusters', 'N/A')}\n\n"

    msg += "工具评估结果:\n"
    for tool, result in quality_results.items():
        msg += f"- {tool}: {result.get('details', 'N/A')} → {result.get('recommendation', 'N/A')}\n"

    msg += "\n请选择:\n"
    msg += "1. 接受当前聚类，继续下一步\n"
    msg += "2. 提高 resolution 重新聚类\n"
    msg += "3. 降低 resolution 重新聚类\n"
    msg += "4. 手动指定 resolution\n"

    return msg
