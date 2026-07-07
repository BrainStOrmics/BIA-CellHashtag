# Configuration

[Deep Agents Code](/oss/python/deepagents/code/overview)

# Configuration
Copy page

Configure Deep Agents Code with config.toml, hooks, and MCP serversCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Deep Agents Code stores its configuration in the `~/.deepagents/` directory. The main config files are:

| 
 | File | Format | Purpose
 | `config.toml` | TOML | Model defaults, provider settings, constructor params, profile overrides, themes, update settings, MCP trust store
 | `.env` | Dotenv | Global API keys and secrets
 | `hooks.json` | JSON | External tool subscriptions to Deep Agents Code lifecycle events
 | `.mcp.json` | JSON | Global MCP server definitions
Files under `~/.deepagents/.state/` hold per-machine Deep Agents Code state and are managed automatically.

## [​](#provider-credentials)Provider credentials

For every provider, Deep Agents Code checks two credential sources in order:

{' ' * (self.list_depth - 1)}- **Stored keys** — entered via `/auth` and persisted locally.

{' ' * (self.list_depth - 1)}- **Environment variables** — `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `GOOGLE_API_KEY`, etc., plus the same names with a [`DEEPAGENTS_CODE_`](#deepagents_code_-prefix) prefix and any keys loaded from [`.env` files](#environment-variables).

A stored key always wins over an env-var key for the same provider.

### [​](#use-/auth-recommended)Use `/auth` (recommended)

Open the credential manager from any session:

```
/auth

```

The manager lists every LLM provider Deep Agents Code knows about, marks the ones that already have a key, and lets you paste, rotate, or clear credentials inline. Each provider renders as one of:

| 
 | Label | Meaning
 | `✓ credentials set` | A stored or env-var key is available
 | `! missing OPENAI_API_KEY` | The named env var is unset and no key is stored — pick the row to paste one
 | `local provider` | No credentials needed (e.g., local Ollama)
 | `implicit auth` | Credentials are implicit (e.g., Vertex AI Application Default Credentials)
 | `custom auth` | A `class_path` provider managing its own auth (mTLS, JWT, custom headers)
Press Enter on a row to paste a key, or use the row’s delete affordance to clear an existing one. Saved keys are written atomically and persist across sessions and `/reload`.
Keys are scoped to your user account on this machine — Deep Agents Code never transmits them anywhere except to the configured provider’s API.
`/auth` only manages **LLM provider** credentials. Tool credentials such as `TAVILY_API_KEY` (web search) and `LANGSMITH_API_KEY` (tracing) are read from the environment instead — [set them in `~/.deepagents/.env` or your shell](#environment-variables).

### [​](#environment-variables-ci-and-headless)Environment variables (CI and headless)

For non-interactive runs, CI/CD pipelines, or anywhere a TUI isn’t available, export the provider’s env var or write it to `~/.deepagents/.env`:

```
export ANTHROPIC_API_KEY="sk-ant-..."
export OPENAI_API_KEY="sk-..."

```

See [Environment variables](#environment-variables) below for loading order, the `DEEPAGENTS_CODE_` prefix, and per-project overrides.

## [​](#environment-variables)Environment variables

Deep Agents Code loads environment variables from dotenv files so you don’t need to `export` API keys in your shell profile or duplicate `.env` files across projects.

### [​](#loading-order-and-precedence)Loading order and precedence

Two `.env` files are loaded at startup:

{' ' * (self.list_depth - 1)}- **Project `.env`** — the `.env` file in your current working directory (if present)

{' ' * (self.list_depth - 1)}- **Global `~/.deepagents/.env`** — a single shared file that acts as a fallback for all projects

The effective precedence is: **shell environment > project `.env` > global `.env`**. Values already set in the shell are never overwritten—including on `/reload`.

### [​](#deepagents_code_-prefix)`DEEPAGENTS_CODE_` prefix

All Deep Agents Code-specific environment variables use a `DEEPAGENTS_CODE_` prefix (e.g., `DEEPAGENTS_CODE_AUTO_UPDATE`, `DEEPAGENTS_CODE_DEBUG`). See the [environment variable reference](#environment-variable-reference) for the full list.
The prefix also works as an override mechanism for any environment variable Deep Agents Code reads, including third-party credentials. Deep Agents Code checks `DEEPAGENTS_CODE_{NAME}` first, then falls back to `{NAME}`:
~/.deepagents/.env
```
# Override OPENAI_API_KEY only for Deep Agents Code, without affecting other tools
DEEPAGENTS_CODE_OPENAI_API_KEY=sk-cli-only

# Block a shell-exported key within Deep Agents Code by setting the prefixed var to empty
DEEPAGENTS_CODE_ANTHROPIC_API_KEY=

```

On `/reload`, Deep Agents Code re-reads `.env` files and picks up prefixed values, so you can rotate keys without restarting.

### [​](#example)Example

Store API keys once in `~/.deepagents/.env`:

```
ANTHROPIC_API_KEY=sk-ant-...
OPENAI_API_KEY=sk-...
LANGSMITH_API_KEY=lsv2_...
TAVILY_API_KEY=tvly-...

# If OPENAI_API_KEY is already exported in your shell for other tools,
# use the prefix to give Deep Agents Code its own key without conflict
DEEPAGENTS_CODE_OPENAI_API_KEY=sk-cli-only-...

```

Then override per-project where needed by placing a `.env` in the project directory.

### [​](#enable-web-search-with-tavily)Enable web search with Tavily

The built-in `web_search` tool uses [Tavily](https://tavily.com). Deep Agents Code shows a “Web search disabled — `TAVILY_API_KEY` is not set” notification on startup until you provide a key. Unlike model provider keys, Tavily is **not** managed by `/auth`; it is read directly from the environment.
1[](#)

Get a keySign up at [tavily.com](https://tavily.com) and copy the key (it starts with `tvly-`). The free tier is sufficient for most Deep Agents Code usage.2[](#)

Add it to your environmentAdd the key to `~/.deepagents/.env` so every session picks it up:~/.deepagents/.env
```
TAVILY_API_KEY=tvly-...

```
Or export it in your shell for a one-off session:
```
export TAVILY_API_KEY=tvly-...

```
Shell exports take precedence over `.env` values (see [Loading order and precedence](#loading-order-and-precedence)). To scope a key to Deep Agents Code only without affecting other tools that read `TAVILY_API_KEY`, use the [`DEEPAGENTS_CODE_` prefix](#deepagents_code_-prefix): `DEEPAGENTS_CODE_TAVILY_API_KEY=tvly-...`.3[](#)

Reload or restartIn an existing session, run `/reload` to re-read `.env` files. On the next launch, the “Web search disabled” notification goes away and the agent can call `web_search`.

---

## [​](#config-file)Config file

`~/.deepagents/config.toml` lets you customize model providers, set defaults, and pass extra parameters to model constructors.

### [​](#default-and-recent-model)Default and recent model

```
[models]
default = "ollama:qwen3:4b" # your intentional long-term preference
recent = "google_genai:gemini-3.5-flash" # last /model switch (written automatically)

```

`[models].default` always takes priority over `[models].recent`. The `/model` command only writes to `[models].recent`, so your configured default is never overwritten by mid-session switches. To remove the default, use `/model --default --clear` or delete the `default` key from the config file.

### [​](#default-and-recent-agent)Default and recent agent

```
[agents]
default = "backend-dev" # your intentional long-term preference (Ctrl+S in /agents picker)
recent = "frontend-dev" # last /agents switch (written automatically)

```

`[agents].default` always takes priority over `[agents].recent`. Selecting an agent in the `/agents` picker with `Enter` writes to `recent`; pressing `Ctrl+S` on the highlighted row pins it as `default`. Pressing `Ctrl+S` again on the same row clears the default.
Explicit `-a`/`--agent` always overrides both, and `-r`/`--resume` bypasses both so the thread’s original agent is restored. See [Command reference](/oss/python/deepagents/code/overview#command-reference) for related flags.

### [​](#provider-configuration)Provider configuration

Each provider is a TOML table under `[models.providers]`:

```
[models.providers.<name>]
models = ["gpt-4o"]
api_key_env = "OPENAI_API_KEY"
base_url = "https://api.openai.com/v1"
class_path = "my_package.models:MyChatModel"
enabled = true

[models.providers.<name>.params]
temperature = 0
max_tokens = 4096

[models.providers.<name>.params."gpt-4o"]
temperature = 0.7

```

Providers have the following configuration options:
[​](#param-models)modelsstring[]optionalA list of model names to show in the interactive `/model` switcher for this provider. For providers that already ship with model profiles, any names you add here appear alongside the bundled ones, useful for newly released models that haven’t been added to the package yet. For [arbitrary providers](#arbitrary-providers), this list is the only source of models in the switcher.Models listed here **bypass** the profile-based [filtering criteria](/oss/python/deepagents/code/providers#which-models-appear-in-the-switcher) and always appear in the switcher. This makes it the recommended way to surface models that are excluded because their profile lacks `tool_calling` support or doesn’t exist yet.This key is optional. You can always pass any model name directly to `/model` or `--model` regardless of whether it appears in the switcher; the provider validates the name at request time.
[​](#param-api-key-env)api_key_envstringoptionalThe **name** of the environment variable that holds the API key (e.g., `"OPENAI_API_KEY"`), not the key itself. Deep Agents Code reads the credential from this env var at startup to verify access before creating the model. Most chat model packages read from a default env var automatically. See the [Provider reference](/oss/python/deepagents/code/providers#provider-reference) table for which variable each provider checks.
[​](#param-base-url)base_urlstringoptionalOverride the base URL used by the provider, if supported. Refer to your provider packages’ [reference docs](https://reference.langchain.com/python/integrations/) for more info.
[​](#param-params)paramsobjectoptionalExtra keyword arguments forwarded to the model constructor. Flat keys (e.g., `temperature = 0`) apply to every model from this provider. Model-keyed sub-tables (e.g., `[params."gpt-4o"]`) override individual values for that model only; the merge is shallow (model wins on conflict).Do not put credentials (e.g., `api_key`) in `params`. Use [`api_key_env`](#provider-configuration) to point at an environment variable instead.
[​](#param-profile)profileobjectoptional(Advanced) Override fields in the model’s runtime [profile](/oss/python/langchain/models#model-profiles) (e.g., `max_input_tokens`). Flat keys apply to every model from this provider. Model-keyed sub-tables (e.g., `[profile."claude-sonnet-4-5"]`) override individual values for that model only; the merge is shallow (model wins on conflict). These overrides are applied after the model is created, so they take effect for context-limit display, auto-summarization, and any other feature that reads the profile.
[​](#param-class-path)class_pathstringoptionalUsed for [arbitrary model](#arbitrary-providers) providers. Fully-qualified Python class in `module.path:ClassName` format. When set, Deep Agents Code imports and instantiates this class directly for provider `<name>`. The class must be a `BaseChatModel` subclass.
[​](#param-enabled)enabledbooleandefault:"true"optionalWhether this provider appears in the `/model` selector. Set to `false` to hide a provider that was auto-discovered from an installed package (e.g., a transitive dependency you don’t want cluttering the switcher). You can still use a disabled provider directly via `/model provider:model` or `--model`.

### [​](#model-constructor-params)Model constructor params

Any provider can use the `params` table to pass extra arguments to the model constructor:

```
[models.providers.ollama.params]
temperature = 0
num_ctx = 8192

```

#### [​](#per-model-overrides)Per-model overrides

If a specific model needs different params, add a model-keyed sub-table under `params` to override individual values without duplicating the entire provider config:

```
[models.providers.ollama]
models = ["qwen3:4b", "llama3"]

[models.providers.ollama.params]
temperature = 0
num_ctx = 8192

[models.providers.ollama.params."qwen3:4b"]
temperature = 0.5
num_ctx = 4000

```

With this configuration:

{' ' * (self.list_depth - 1)}- `ollama:qwen3:4b` gets `{temperature: 0.5, num_ctx: 4000}` — model overrides win.

{' ' * (self.list_depth - 1)}- `ollama:llama3` gets `{temperature: 0, num_ctx: 8192}` — no override, provider-level params only.

The merge is shallow: any key present in the model sub-table replaces the same key from the provider-level params, while keys only at the provider level are preserved.
For one-off adjustments without editing `config.toml`, pass a JSON object via `--model-params` at launch or mid-session with `/model`. CLI flags take highest priority over the config file. See [Model parameters](/oss/python/deepagents/code/providers#model-parameters) on the providers page for syntax and provider-specific examples.

### [​](#profile-overrides-advanced)Profile overrides (Advanced)

Override fields in the model’s runtime profile to change how Deep Agents Code interprets model capabilities. The most common use case is lowering `max_input_tokens` to trigger auto-summarization earlier — useful for testing or for constraining context usage:

```
# Apply to all models from this provider
[models.providers.anthropic.profile]
max_input_tokens = 4096

```

Per-model sub-tables work the same way as `params` — the model-level value wins on conflict:

```
[models.providers.anthropic.profile]
max_input_tokens = 4096

# This model gets a higher limit
[models.providers.anthropic.profile."claude-sonnet-4-5"]
max_input_tokens = 8192

```

Profile overrides are merged into the model’s profile after creation. Any feature that reads the profile — context-limit display in the status bar, auto-summarization thresholds, capability checks — will see the overridden values.

#### [​](#cli-profile-overrides-with-profile-override-advanced)CLI profile overrides with `--profile-override` (Advanced)

To override model profile fields at runtime without editing the config file, pass a JSON object via `--profile-override`:

```
dcode --profile-override '{"max_input_tokens": 4096}'

# Combine with --model
dcode --model google_genai:gemini-3.5-flash --profile-override '{"max_input_tokens": 4096}'

# In non-interactive mode
dcode -n "Summarize this repo" --profile-override '{"max_input_tokens": 4096}'

```

These are merged on top of config file profile overrides (CLI wins). The priority chain is: model default < config.toml profile < CLI `--profile-override`.
`--profile-override` values persist across mid-session `/model` hot-swaps — switching models re-applies the override to the new model.

### [​](#custom-base-url)Custom base URL

Some provider packages accept a `base_url` to override the default endpoint. For example, `langchain-ollama` defaults to `http://localhost:11434` via the underlying `ollama` client. To point it elsewhere, set `base_url` in your configuration:

```
[models.providers.ollama]
base_url = "http://your-host-here:port"

```

Refer to your provider’s reference documentation for compatibility information and additional considerations.

### [​](#compatible-apis)Compatible APIs

For providers that expose APIs that are wire-compatible with OpenAI or Anthropic, you can use the existing `langchain-openai` or `langchain-anthropic` packages by pointing `base_url` at the provider’s endpoint:

```
[models.providers.openai]
base_url = "https://api.example.com/v1"
api_key_env = "EXAMPLE_API_KEY"
models = ["my-model"]

```

```
[models.providers.anthropic]
base_url = "https://api.example.com"
api_key_env = "EXAMPLE_API_KEY"
models = ["my-model"]

```

Any features added on top of the official spec by the provider will not be captured. If the provider offers a dedicated LangChain integration package, prefer that instead.
The OpenAI provider defaults to the [Responses API](https://platform.openai.com/docs/api-reference/responses), which most OpenAI-compatible gateways do not implement. If your provider only supports the Chat Completions API, invocation will likely fail. Disable the Responses API explicitly:
```
[models.providers.openai.params]
use_responses_api = false

```

### [​](#adding-models-to-the-interactive-switcher)Adding models to the interactive switcher

Some providers (e.g. `langchain-ollama`) don’t bundle model profile data (see [Provider reference](/oss/python/deepagents/code/providers#provider-reference) for full listing). When this is the case, the interactive `/model` switcher won’t list models for that provider. You can fill in the gap by defining a `models` list in your config file for the provider:

```
[models.providers.ollama]
models = ["llama3", "mistral", "codellama"]

```

The `/model` switcher will now include an Ollama section with these models listed.
This is entirely optional. You can always switch to any model by specifying its full name directly:

```
/model ollama:llama3

```

When `langchain-ollama` is installed and the daemon is reachable, Deep Agents Code auto-discovers locally pulled models and merges them into the switcher—no `models` list required. Run `/reload` to refresh after pulling new models, or set `DEEPAGENTS_CODE_OLLAMA_DISCOVERY=0` to opt out.

### [​](#arbitrary-providers)Arbitrary providers

You can use any [LangChain `BaseChatModel`](https://reference.langchain.com/python/langchain_core/language_models/#langchain_core.language_models.BaseChatModel) subclass using `class_path`. Deep Agents Code imports and instantiates the class directly — no built-in provider package required.

```
[models.providers.my_custom]
class_path = "my_package.models:MyChatModel"
api_key_env = "MY_API_KEY"
base_url = "https://my-endpoint.example.com"

[models.providers.my_custom.params]
temperature = 0
max_tokens = 4096

```

`api_key_env` and `base_url` are optional. `class_path` providers are expected to handle their own authentication internally — useful when your model uses custom auth (JWT tokens, proprietary headers, mTLS, etc.) rather than a standard API key:

```
[models.providers.xyz]
class_path = "abc.integrations.deepagents:DeepAgentsXYZChat"
models = ["abc-xyz-1"]

[models.providers.xyz.params]
bypass_auth = true
temperature = 0

```

With this config, switch to the model with `/model xyz:abc-xyz-1` or `--model xyz:abc-xyz-1`.
Deep Agents Code requires **tool calling** support. If your custom model supports tool calling but Deep Agents Code doesn’t know about it, declare it in the provider profile:
```
[models.providers.xyz.profile]
tool_calling = true
max_input_tokens = 128000

```
Set `max_input_tokens` to what your model supports to enable accurate context length tracking and auto-summarization.
The provider package must be installed in the same Python environment as `deepagents-code`:

```
# If deepagents-code was installed with uv tool:
uv tool install deepagents-code --with my_package

```

When you switch to `my_custom:my-model-v1` (via `/model` or `--model`), the model name (`my-model-v1`) is passed as the `model` kwarg:

```
MyChatModel(model="my-model-v1", base_url="...", api_key="...", temperature=0, max_tokens=4096)

```

`class_path` executes arbitrary Python code from your config file. This has the same trust model as `pyproject.toml` build scripts — you control your own machine.
Your provider package may optionally provide model profiles at a `_PROFILES` dict in `<package>.data._profiles` in lieu of defining them under the `models` key. See LangChain [model profiles](https://github.com/langchain-ai/langchain/tree/master/libs/model-profiles) for more info.

---

## [​](#skills-extra-allowed-directories)Skills extra allowed directories

By default, when Deep Agents Code loads skills it validates that a resolved skill file path stays inside one of the standard [skill directories](/oss/python/deepagents/code/data-locations#skills). This prevents symlinks inside skill directories from reading arbitrary files outside those roots.
If you store shared skill assets in a non-standard location and use symlinks from a standard skill directory to reference them, you can add that location to the containment allowlist. This does **not** add a new skill discovery location: skills are still only discovered from the standard directories.
[​](#param-extra-allowed-dirs)extra_allowed_dirsstring[]optionalPaths added to the skill containment allowlist. Supports `~` expansion.
```
[skills]
extra_allowed_dirs = [
 "~/shared-skills",
 "/opt/team-skills",
]

```

Alternatively, set the `DEEPAGENTS_CODE_EXTRA_SKILLS_DIRS` environment variable as a colon-separated list:

```
export DEEPAGENTS_CODE_EXTRA_SKILLS_DIRS="~/shared-skills:/opt/team-skills"

```

When the environment variable is set, it takes precedence over the config file value. Changes take effect on `/reload`.

---

## [​](#themes)Themes

Use `/theme` to open an interactive theme selector. Navigate the list to preview themes in real-time, press `Enter` to persist your choice to `config.toml`.
Deep Agents Code ships with many built-in themes. The default theme is `langchain`, a dark theme with LangChain-branded colors. The selected theme is persisted under `[ui]`:

```
[ui]
theme = "langchain-dark"

```

### [​](#user-defined-themes)User-defined themes

Define custom themes under `[themes.<name>]` sections in `config.toml`. Each section requires `label` (str). `dark` (bool) defaults to `false` if omitted — set to `true` for dark themes. All color fields are optional — omitted fields fall back to the built-in dark or light palette based on the `dark` flag.

```
[themes.my-solarized]
label = "My Solarized"
dark = true
primary = "#268BD2"
warning = "#B58900"

# Theme names with spaces require TOML quoting
[themes."ocean breeze"]
label = "Ocean Breeze"
primary = "#0077B6"
background = "#CAF0F8"

```

User-defined themes appear alongside built-in themes in the `/theme` selector.

### [​](#override-built-in-theme-colors)Override built-in theme colors

To tweak a built-in theme’s colors without creating a new theme, use a `[themes.<builtin-name>]` section. Only color fields are read — `label` and `dark` are inherited from the built-in:

```
[themes.langchain]
primary = "#FF5500"

```

Omitted color fields retain the existing built-in values.
Changes to `[themes.*]` sections take effect on `/reload`.

### [​](#map-themes-to-terminals)Map themes to terminals

If you switch between terminals with different color schemes (for example, a dark iTerm and a light Apple Terminal), map each one to a theme under `[ui.terminal_themes]`. Deep Agents Code matches the shell’s `TERM_PROGRAM` and applies the mapped theme automatically:

```
[ui.terminal_themes]
"Apple_Terminal" = "langchain-light"
"iTerm.app" = "langchain"

```

Press `T` in the `/theme` picker to save the highlighted theme for the current terminal, or run `echo $TERM_PROGRAM` to find your terminal’s identifier and add it by hand.

Advanced: picker shortcuts, resolution order, terminal identifiers

#### [​](#picker-shortcuts)Picker shortcuts
In the `/theme` selector:

{' ' * (self.list_depth - 1)}- `N` toggles between display labels and canonical registry keys—the keys are what `[ui] theme` and `[ui.terminal_themes]` accept.

{' ' * (self.list_depth - 1)}- `T` saves the highlighted theme into `[ui.terminal_themes]` for the current `TERM_PROGRAM`. The mapped theme is badged `(default)` in the picker.

#### [​](#common-term_program-values)Common `TERM_PROGRAM` values
Keys are matched verbatim against the environment variable—quote them in TOML when they contain dots or special characters.

| 
 | Terminal | `TERM_PROGRAM`
 | Apple Terminal | `Apple_Terminal`
 | iTerm2 | `iTerm.app`
 | WezTerm | `WezTerm`
 | VS Code integrated terminal | `vscode`
 | Ghostty | `ghostty`

#### [​](#resolution-order)Resolution order
Deep Agents Code resolves a theme on every launch using this precedence:

{' ' * (self.list_depth - 1)}- `DEEPAGENTS_CODE_THEME` environment variable (explicit override).

{' ' * (self.list_depth - 1)}- `[ui.terminal_themes]` mapping for the current `TERM_PROGRAM`.

{' ' * (self.list_depth - 1)}- `[ui] theme` saved preference (set by `/theme`).

{' ' * (self.list_depth - 1)}- The built-in default (`langchain`).

---

## [​](#auto-update)Auto-update

Deep Agents Code can automatically check for and install updates.

{' ' * (self.list_depth - 1)}- Config file
{' ' * (self.list_depth - 1)}- Environment variable
```
[update]
auto_update = true

```

```
export DEEPAGENTS_CODE_AUTO_UPDATE=1

```

The environment variable takes precedence over the config file.
When enabled, Deep Agents Code checks PyPI for a newer version at session start and automatically upgrades using the detected install method (uv, Homebrew, or pip). When disabled (default), Deep Agents Code shows an update hint with the appropriate install command instead.
You can also check for and install updates manually at any time with the `/update` slash command, which bypasses the cache and reports success or failure inline.
After an upgrade, Deep Agents Code shows a “what’s new” banner on the next launch with a link to the changelog.
At session exit, if a newer version was detected during the session, an update banner is displayed as a reminder.

---

## [​](#managed-deployments)Managed deployments

The [install script](https://github.com/langchain-ai/deepagents/blob/main/libs/cli/scripts/install.sh) supports running as root, targeting macOS MDM tools (Kandji, Jamf, etc.) that execute scripts in a minimal root environment.
When `id -u` is `0`, the script:

{' ' * (self.list_depth - 1)}- Resolves the real console user’s `HOME` (via `/dev/console` or a `/Users` directory scan)

{' ' * (self.list_depth - 1)}- `chown`s all created files back to the target user after each install step

Non-root installs are unaffected: all root-specific code paths short-circuit when not running as root.
To pre-configure auto-update for managed installs, set `DEEPAGENTS_CODE_AUTO_UPDATE=1` in the user’s shell profile or deploy a `config.toml` with `[update] auto_update = true` to `~/.deepagents/config.toml`. To suppress automatic updates and update checks entirely, set `DEEPAGENTS_CODE_NO_UPDATE_CHECK=1`.

---

## [​](#environment-variable-reference)Environment variable reference

All Deep Agents Code-specific environment variables use the `DEEPAGENTS_CODE_` prefix. See [`DEEPAGENTS_CODE_` prefix](#deepagents_code_-prefix) for how the prefix also works as an override for third-party credentials.
[​](#param-deepagents-code-auto-update)DEEPAGENTS_CODE_AUTO_UPDATEstringoptionalEnable automatic Deep Agents Code updates (`1`, `true`, or `yes`).
[​](#param-deepagents-code-debug)DEEPAGENTS_CODE_DEBUGstringoptionalEnable verbose debug logging to a file. Accepts `1`, `true`, `yes`, `on` (case-insensitive) as enabled; `0`, `false`, `no`, `off`, empty string, or unset disables it. When enabled, the per-session server log file is preserved on shutdown and its path is printed to stderr for triage.
[​](#param-deepagents-code-debug-file)DEEPAGENTS_CODE_DEBUG_FILEstringdefault:"/tmp/deepagents_debug.log"optionalPath for the debug log file.
[​](#param-deepagents-code-extra-skills-dirs)DEEPAGENTS_CODE_EXTRA_SKILLS_DIRSstringoptionalColon-separated paths added to the [skill containment allowlist](#skills-extra-allowed-directories).
[​](#param-deepagents-code-langsmith-project)DEEPAGENTS_CODE_LANGSMITH_PROJECTstringoptionalOverride the LangSmith project name for agent traces. See [Trace with LangSmith](/oss/python/deepagents/code/overview#trace-with-langsmith).
[​](#param-deepagents-code-no-update-check)DEEPAGENTS_CODE_NO_UPDATE_CHECKstringoptionalDisable automatic update checking when set.
[​](#param-deepagents-code-shell-allow-list)DEEPAGENTS_CODE_SHELL_ALLOW_LISTstringoptionalComma-separated shell commands to allow (or `recommended` / `all`).
[​](#param-deepagents-code-user-id)DEEPAGENTS_CODE_USER_IDstringoptionalAttach a user identifier to LangSmith trace metadata.

---

## [​](#external-editor)External editor

Press `Ctrl+X` or type `/editor` to compose prompts in an external editor. Deep Agents Code checks `$VISUAL`, then `$EDITOR`, then falls back to `vi` (macOS/Linux) or `notepad` (Windows). GUI editors (VS Code, Cursor, Zed, Sublime Text, Windsurf) automatically receive a `--wait` flag so Deep Agents Code blocks until you close the file.

```
# Set in your shell profile (~/.zshrc, ~/.bashrc, etc.)
export VISUAL="code" # GUI editor (--wait auto-injected)
export EDITOR="nvim" # Terminal fallback

```

---

## [​](#hooks)Hooks

Hooks let external programs react to Deep Agents Code lifecycle events. Configure commands in `~/.deepagents/hooks.json` and it pipes a JSON payload to each matching command’s stdin whenever an event fires.
Hooks run fire-and-forget in a background thread — they never block Deep Agents Code and failures are logged without interrupting your session.

### [​](#setup)Setup

Create `~/.deepagents/hooks.json`:

```
{
 "hooks": [
 {
 "command": ["bash", "-c", "cat >> ~/deepagents-events.log"],
 "events": ["session.start", "session.end"]
 }
 ]
}

```

Now every time a session starts or ends, Deep Agents Code appends the event payload to `~/deepagents-events.log`.

### [​](#hook-configuration)Hook configuration

The config file contains a single `hooks` array. Each entry has:
[​](#param-command)commandlist[str]requiredCommand and arguments to run. No shell expansion: use `["bash", "-c", "..."]` if needed.
[​](#param-events)eventslist[str]optionalEvent names to subscribe to. Omit or leave empty to receive **all** events.

```
{
 "hooks": [
 {
 "command": ["python3", "my_handler.py"],
 "events": ["session.start", "task.complete"]
 },
 {
 "command": ["bash", "log_everything.sh"]
 }
 ]
}

```

The second hook above has no `events` filter, so it receives every event Deep Agents Code emits.

### [​](#payload-format)Payload format

Each hook command receives a JSON object on stdin with an `"event"` key plus event-specific fields:

```
{
 "event": "session.start",
 "thread_id": "abc123"
}

```

### [​](#events-reference)Events reference

#### [​](#session-start)`session.start`

Fired when an agent session begins (both interactive and non-interactive modes).
[​](#param-thread-id)thread_idstringrequiredThe session thread identifier.

#### [​](#session-end)`session.end`

Fired when a session exits.
[​](#param-thread-id-1)thread_idstringrequiredThe session thread identifier.

#### [​](#user-prompt)`user.prompt`

Fired in interactive mode when the user submits a chat message.
No additional fields.

#### [​](#input-required)`input.required`

Fired when the agent requires human input (human-in-the-loop interrupt).
No additional fields.

#### [​](#permission-request)`permission.request`

Fired before the approval dialog when one or more tool calls need user permission.
[​](#param-tool-names)tool_nameslist[str]requiredNames of the tools requesting approval.

#### [​](#tool-error)`tool.error`

Fired when a tool call returns an error.
[​](#param-tool-names-1)tool_nameslist[str]requiredNames of the tool(s) that errored.

#### [​](#task-complete)`task.complete`

Fired when the agent finishes its current task (the streaming loop ends without further interrupts).
[​](#param-thread-id-2)thread_idstringrequiredThe session thread identifier.

#### [​](#context-compact)`context.compact`

Fired before Deep Agents Code compacts (summarizes) the conversation context.
No additional fields.

### [​](#execution-model)Execution model

{' ' * (self.list_depth - 1)}- **Background thread**: Hook subprocesses run in a thread via `asyncio.to_thread` so the main event loop is never blocked.

{' ' * (self.list_depth - 1)}- **Concurrent dispatch**: When multiple hooks match an event, they run concurrently in a thread pool.

{' ' * (self.list_depth - 1)}- **5-second timeout**: Each command has a 5-second timeout. Commands that exceed this are killed.

{' ' * (self.list_depth - 1)}- **Fire-and-forget**: Errors are caught per-hook and logged at debug/warning level. A failing hook never crashes or stalls Deep Agents Code.

{' ' * (self.list_depth - 1)}- **Lazy loading**: The config file is read once on the first event dispatch and cached for the rest of the session.

{' ' * (self.list_depth - 1)}- **No shell expansion**: Commands are executed directly (not through a shell). Wrap in `["bash", "-c", "..."]` if you need shell features like pipes or variable expansion.

### [​](#hook-examples)Hook examples

Log all events to a file
```
{
 "hooks": [
 {
 "command": ["bash", "-c", "jq -c . >> ~/.deepagents/hook-events.jsonl"],
 "events": []
 }
 ]
}

```

Desktop notification on task completion (macOS)
```
{
 "hooks": [
 {
 "command": [
 "bash", "-c",
 "osascript -e 'display notification \"Agent finished\" with title \"Deep Agents\"'"
 ],
 "events": ["task.complete"]
 }
 ]
}

```

Python handlerWrite a handler script that reads the JSON payload from stdin:my_handler.py
```
import json
import sys

payload = json.load(sys.stdin)
event = payload["event"]

if event == "session.start":
 print(f"Session started: {payload['thread_id']}", file=sys.stderr)
elif event == "permission.request":
 print(f"Approval needed for: {payload['tool_names']}", file=sys.stderr)

```
~/.deepagents/hooks.json
```
{
 "hooks": [
 {
 "command": ["python3", "my_handler.py"],
 "events": ["session.start", "permission.request"]
 }
 ]
}

```

### [​](#security-considerations)Security considerations

Hooks follow the same trust model as Git hooks or shell aliases — any user who can write to `~/.deepagents/hooks.json` can execute arbitrary commands. This is by design:

{' ' * (self.list_depth - 1)}- **No command injection**: Payload data flows only to stdin as JSON, never to command-line arguments. `json.dumps` handles escaping.

{' ' * (self.list_depth - 1)}- **No shell by default**: Commands run with `shell=False`, preventing shell injection.

{' ' * (self.list_depth - 1)}- **Malformed config**: Invalid JSON or unexpected types produce logged warnings, not security issues.

Only add hooks from sources you trust. A hook has the same permissions as your user account.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/configuration.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Model providersPrevious](/oss/python/deepagents/code/providers)[MCP toolsNext](/oss/python/deepagents/code/mcp-tools)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
